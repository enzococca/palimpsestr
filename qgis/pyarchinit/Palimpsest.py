# -*- coding: utf-8 -*-
"""
Palimpsest Analysis Dialog for pyArchInit.

Runs the palimpsestr Stratigraphic Entanglement Field (SEF) analysis directly
on the currently connected pyArchInit database, via the QGIS Processing R
Provider algorithms shipped in the palimpsestr `qgis/processing` folder
(`r:fit_sef_from_pyarchinit_db`, `r:detect_intrusions_from_pyarchinit_db`).

Modelled on tabs/Movecost.py. Drop this file into the pyArchInit `tabs/` folder
and wire it from pyarchinitPlugin.py (see the snippet at the bottom of this
file).

@author: Enzo Cocca <enzo.ccc@gmail.com>
"""
import os
from urllib.parse import urlparse

from qgis.PyQt.QtWidgets import (
    QDialog, QVBoxLayout, QFormLayout, QHBoxLayout, QLabel, QSpinBox,
    QDoubleSpinBox, QComboBox, QCheckBox, QLineEdit, QPushButton, QMessageBox)
from qgis.core import (QgsApplication, QgsCategorizedSymbolRenderer,
                       QgsRendererCategory, QgsSymbol)
from qgis.PyQt.QtGui import QColor
import processing

# Okabe-Ito colour-blind-safe palette, matching palimpsestr's .phase_colours().
PHASE_COLOURS = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7",
                 "#D55E00", "#0072B2", "#F0E442", "#999999"]

FIT_ALG = "r:fit_sef_from_pyarchinit_db"
INTRUSIONS_ALG = "r:detect_intrusions_from_pyarchinit_db"


class pyarchinit_Palimpsest(QDialog):
    HOME = os.environ.get('PYARCHINIT_HOME', os.path.expanduser("~"))

    def __init__(self, iface):
        super().__init__()
        self.iface = iface
        self.setWindowTitle("palimpsestr — Palimpsest analysis")
        self._build_ui()

    # ------------------------------------------------------------------ UI ---
    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "Probabilistic palimpsest decomposition (SEF) on the active "
            "pyArchInit database."))

        form = QFormLayout()
        self.spin_k = QSpinBox(); self.spin_k.setRange(2, 12); self.spin_k.setValue(4)
        form.addRow("Number of phases (K):", self.spin_k)

        self.combo_model = QComboBox()
        self.combo_model.addItems(["multinomial (recommended)", "gaussian (legacy)"])
        form.addRow("Class model:", self.combo_model)

        self.check_noise = QCheckBox("Noise / outlier component"); self.check_noise.setChecked(True)
        form.addRow("", self.check_noise)

        self.spin_thr = QDoubleSpinBox(); self.spin_thr.setRange(0.0, 1.0)
        self.spin_thr.setSingleStep(0.05); self.spin_thr.setValue(0.5)
        form.addRow("Intrusion threshold:", self.spin_thr)

        self.edit_site = QLineEdit(); self.edit_site.setPlaceholderText("(all sites)")
        form.addRow("Site filter (optional):", self.edit_site)
        layout.addLayout(form)

        self.lbl_db = QLabel(self._describe_db())
        self.lbl_db.setWordWrap(True)
        layout.addWidget(self.lbl_db)

        buttons = QHBoxLayout()
        self.btn_fit = QPushButton("Fit SEF model"); self.btn_fit.clicked.connect(self.run_fit)
        self.btn_intr = QPushButton("Detect intrusions"); self.btn_intr.clicked.connect(self.run_intrusions)
        buttons.addWidget(self.btn_fit); buttons.addWidget(self.btn_intr)
        layout.addLayout(buttons)

    # ----------------------------------------------------- active DB params ---
    def _active_conn_str(self):
        """Return the SQLAlchemy connection string of the active pyArchInit DB."""
        try:
            from ..modules.db.pyarchinit_conn_strings import Connection
        except Exception:
            try:
                from modules.db.pyarchinit_conn_strings import Connection
            except Exception:
                return None
        try:
            return Connection().conn_str()
        except Exception:
            return None

    def _db_params(self):
        """Build the connection parameters expected by the .rsx algorithms."""
        params = {'Database_file': '', 'PostgreSQL_host': '', 'PostgreSQL_dbname': '',
                  'PostgreSQL_user': '', 'PostgreSQL_password': ''}
        cs = self._active_conn_str()
        if not cs:
            return params
        if cs.startswith('sqlite'):
            params['Database_file'] = cs.split('sqlite:///', 1)[-1]
        else:
            u = urlparse(cs)
            params['PostgreSQL_host'] = u.hostname or ''
            params['PostgreSQL_dbname'] = (u.path or '').lstrip('/')
            params['PostgreSQL_user'] = u.username or ''
            params['PostgreSQL_password'] = u.password or ''
        return params

    def _describe_db(self):
        p = self._db_params()
        if p['Database_file']:
            return "Active database (SQLite/Spatialite): %s" % p['Database_file']
        if p['PostgreSQL_host']:
            return "Active database (PostgreSQL): %s@%s/%s" % (
                p['PostgreSQL_user'], p['PostgreSQL_host'], p['PostgreSQL_dbname'])
        return "No active pyArchInit database connection detected."

    # ----------------------------------------------------------- run algos ---
    def _check_provider(self):
        if QgsApplication.processingRegistry().algorithmById(FIT_ALG) is None:
            QMessageBox.warning(
                self, "palimpsestr",
                "The palimpsestr Processing R scripts are not registered.\n\n"
                "Install the 'Processing R Provider' plugin and copy the "
                "palimpsestr qgis/processing/*.rsx files into the Processing R "
                "scripts folder, then restart QGIS.")
            return False
        return True

    def _site(self):
        return self.edit_site.text().strip()

    def run_fit(self):
        if not self._check_provider():
            return
        params = self._db_params()
        params.update({
            'Site': self._site(),
            'K': self.spin_k.value(),
            'Class_model': self.combo_model.currentIndex(),   # 0 multinomial, 1 gaussian
            'Noise': self.check_noise.isChecked(),
            'Phases': 'TEMPORARY_OUTPUT',
            'Links': 'TEMPORARY_OUTPUT',
            'Diagnostics': 'TEMPORARY_OUTPUT'})
        try:
            res = processing.runAndLoadResults(FIT_ALG, params)
        except Exception as e:
            QMessageBox.critical(self, "palimpsestr", "Analysis failed:\n%s" % e)
            return
        self._style_phases(res.get('Phases'))
        QMessageBox.information(self, "palimpsestr", "SEF analysis complete. Layers loaded.")

    def run_intrusions(self):
        if not self._check_provider():
            return
        params = self._db_params()
        params.update({
            'Site': self._site(),
            'K': self.spin_k.value(),
            'Threshold': self.spin_thr.value(),
            'Intrusions': 'TEMPORARY_OUTPUT'})
        try:
            processing.runAndLoadResults(INTRUSIONS_ALG, params)
        except Exception as e:
            QMessageBox.critical(self, "palimpsestr", "Analysis failed:\n%s" % e)
            return
        QMessageBox.information(self, "palimpsestr", "Intrusion detection complete. Layer loaded.")

    # ----------------------------------------------------------- styling ----
    def _style_phases(self, layer_id):
        from qgis.core import QgsProject
        layer = None
        if layer_id:
            found = QgsProject.instance().mapLayer(layer_id) if isinstance(layer_id, str) else layer_id
            layer = found
        if layer is None:
            return
        try:
            field = "dominant_phase"
            phases = sorted({f[field] for f in layer.getFeatures() if f[field] is not None})
            categories = []
            for i, ph in enumerate(phases):
                sym = QgsSymbol.defaultSymbol(layer.geometryType())
                sym.setColor(QColor(PHASE_COLOURS[i % len(PHASE_COLOURS)]))
                categories.append(QgsRendererCategory(ph, sym, "phase %s" % ph))
            layer.setRenderer(QgsCategorizedSymbolRenderer(field, categories))
            layer.triggerRepaint()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Wiring in pyarchinitPlugin.py (mirroring actionMovecost), e.g. in initGui():
#
#   from .tabs.Palimpsest import pyarchinit_Palimpsest
#   icon = QIcon(':/plugins/pyarchinit/resources/...')   # any icon
#   self.actionPalimpsest = QAction(icon, "palimpsestr - Analisi palinsesti",
#                                    self.iface.mainWindow())
#   self.actionPalimpsest.triggered.connect(self.runPalimpsest)
#   self.analysisToolButton.addActions([self.actionPalimpsest])
#
#   def runPalimpsest(self):
#       dlg = pyarchinit_Palimpsest(self.iface)
#       dlg.exec_()
# ---------------------------------------------------------------------------
