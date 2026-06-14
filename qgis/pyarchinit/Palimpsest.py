# -*- coding: utf-8 -*-
"""
Palimpsest Analysis Dialog for pyArchInit.

Runs the palimpsestr Stratigraphic Entanglement Field (SEF) analysis directly
on the currently connected pyArchInit SQLite/Spatialite database, via the QGIS
Processing R Provider algorithms (r:palimpsestrfit, r:palimpsestrintrusions),
which the dialog installs/updates into the Processing R scripts folder itself.

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
                       QgsRendererCategory, QgsSymbol, QgsProject)
from qgis.PyQt.QtGui import QColor
import processing

PHASE_COLOURS = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7",
                 "#D55E00", "#0072B2", "#F0E442", "#999999"]

FIT_ALG = "r:palimpsestrfit"
INTRUSIONS_ALG = "r:palimpsestrintrusions"

# Processing R scripts shipped with palimpsestr, embedded so the dialog can
# install/update them itself.
RSX_SCRIPTS = {
    "palimpsestr_fit_db.rsx": r"""##palimpsestr=group
##Palimpsestr Fit=name
##Database_file=file
##Site=string all
##K=number 4
##Class_model=enum literal multinomial;gaussian
##Noise=boolean True
##Phases=output vector
##Links=output vector
##Diagnostics=output table

# Probabilistic palimpsest decomposition straight from a pyArchInit SQLite/
# Spatialite database. Reads inventario_materiali + us (+ US polygon geometry)
# via read_pyarchinit(), fits the Stratigraphic Entanglement Field model, and
# returns a phase-assignment point layer, a high-SEI link layer, and a
# diagnostics table.
library(palimpsestr)
library(sf)
library(DBI)

con  <- DBI::dbConnect(RSQLite::SQLite(), Database_file)
geom <- tryCatch(sf::st_read(Database_file, layer = "pyunitastratigrafiche", quiet = TRUE),
                 error = function(e) NULL)

site <- if (exists("Site") && nchar(Site) > 0 && Site != "all") Site else NULL
d <- read_pyarchinit(con, us_geometry = geom, sito = site)
DBI::dbDisconnect(con)

class_model <- if (is.numeric(Class_model)) c("multinomial", "gaussian")[Class_model + 1] else as.character(Class_model)
use_noise   <- isTRUE(as.logical(Noise))

fit <- fit_sef(d, k = as.integer(K), context = "context",
               tafonomy = "taf_score", class_model = class_model, noise = use_noise)
fit <- reorder_phases(fit)

crs_val <- if (!is.null(geom)) sf::st_crs(geom)$epsg else NA_integer_
if (is.null(crs_val) || is.na(crs_val)) crs_val <- NA_integer_

Phases      <- as_sf_phase(fit, crs = crs_val)
Links       <- as_sf_links(fit, crs = crs_val)
Diagnostics <- as_phase_table(fit)
""",
    "palimpsestr_intrusions_db.rsx": r"""##palimpsestr=group
##Palimpsestr Intrusions=name
##Database_file=file
##Site=string all
##K=number 4
##Threshold=number 0.5
##Intrusions=output vector

# Model-based intrusion detection straight from a pyArchInit SQLite/Spatialite
# database. Fits the SEF model with a noise component and returns the finds as
# a point layer carrying the outlier posterior (intrusion_prob), the
# chronological direction, and the intrusion_type classification.
library(palimpsestr)
library(sf)
library(DBI)

con  <- DBI::dbConnect(RSQLite::SQLite(), Database_file)
geom <- tryCatch(sf::st_read(Database_file, layer = "pyunitastratigrafiche", quiet = TRUE),
                 error = function(e) NULL)

site <- if (exists("Site") && nchar(Site) > 0 && Site != "all") Site else NULL
d <- read_pyarchinit(con, us_geometry = geom, sito = site)
DBI::dbDisconnect(con)

fit <- fit_sef(d, k = as.integer(K), context = "context",
               tafonomy = "taf_score", noise = TRUE)
fit <- reorder_phases(fit)
di <- detect_intrusions(fit, intrusion_threshold = Threshold)

crs_val <- if (!is.null(geom)) sf::st_crs(geom)$epsg else NA_integer_
if (is.null(crs_val) || is.na(crs_val)) crs_val <- NA_integer_

pts <- as_sf_phase(fit, crs = crs_val)
pts$intrusion_prob <- di$intrusion_prob
pts$direction      <- as.character(di$direction)
pts$intrusion_type <- as.character(di$intrusion_type)
Intrusions <- pts
""",
}


class pyarchinit_Palimpsest(QDialog):
    HOME = os.environ.get('PYARCHINIT_HOME', os.path.expanduser("~"))

    def __init__(self, iface):
        super().__init__()
        self.iface = iface
        self.setWindowTitle("palimpsestr \u2014 Palimpsest analysis")
        self._build_ui()
        self.install_scripts(silent=True)

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
        self.btn_install = QPushButton("Install/update R scripts")
        self.btn_install.clicked.connect(lambda: self.install_scripts(silent=False))
        buttons.addWidget(self.btn_fit); buttons.addWidget(self.btn_intr); buttons.addWidget(self.btn_install)
        layout.addLayout(buttons)

    # ----------------------------------------------------- active DB info ---
    def _active_conn_str(self):
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

    def _sqlite_path(self):
        """Return the SQLite/Spatialite DB path, or None if not a SQLite DB."""
        cs = self._active_conn_str()
        if cs and cs.startswith('sqlite'):
            return cs.split('sqlite:///', 1)[-1]
        return None

    def _is_postgres(self):
        cs = self._active_conn_str()
        return bool(cs) and cs.startswith('postgres')

    def _describe_db(self):
        p = self._sqlite_path()
        if p:
            return "Active database (SQLite/Spatialite): %s" % p
        if self._is_postgres():
            return ("Active database is PostgreSQL. These algorithms currently "
                    "read SQLite/Spatialite databases; connect to a SQLite "
                    "pyArchInit database, or run read_pyarchinit() in R.")
        return "No active pyArchInit database connection detected."

    # ----------------------------------------------------------- run algos ---
    def _check_provider(self):
        if QgsApplication.processingRegistry().algorithmById(FIT_ALG) is None:
            QMessageBox.warning(
                self, "palimpsestr",
                "The palimpsestr Processing R scripts are not registered.\n\n"
                "Enable the 'Processing R Provider' plugin and restart QGIS. "
                "The scripts are installed automatically into the Processing R "
                "scripts folder when this dialog opens.")
            return False
        return True

    def _require_sqlite(self):
        path = self._sqlite_path()
        if not path:
            QMessageBox.warning(
                self, "palimpsestr",
                "No active SQLite/Spatialite pyArchInit connection.\n\n"
                "Connect to a SQLite pyArchInit database first.")
            return None
        if not os.path.exists(path):
            QMessageBox.warning(self, "palimpsestr",
                                "Database file not found:\n%s" % path)
            return None
        return path

    def _site(self):
        return self.edit_site.text().strip() or "all"

    def run_fit(self):
        if not self._check_provider():
            return
        path = self._require_sqlite()
        if not path:
            return
        import tempfile
        out = tempfile.mkdtemp(prefix="palimpsestr_")
        ph = os.path.join(out, "sef_phases.gpkg")
        lk = os.path.join(out, "sef_links.gpkg")
        dg = os.path.join(out, "sef_diagnostics.csv")
        params = {
            'Database_file': path,
            'Site': self._site(),
            'K': self.spin_k.value(),
            'Class_model': self.combo_model.currentIndex(),
            'Noise': self.check_noise.isChecked(),
            'Phases': ph, 'Links': lk, 'Diagnostics': dg}
        try:
            res = processing.run(FIT_ALG, params)
        except Exception as e:
            QMessageBox.critical(self, "palimpsestr", "Analysis failed:\n%s" % e)
            return
        n = self._load_outputs([(res.get('Phases', ph), "SEF phases", True),
                                (res.get('Links', lk), "SEF links", False)])
        if n == 0:
            QMessageBox.warning(self, "palimpsestr",
                "The analysis ran but produced no loadable layers (no finds with "
                "coordinates/dating?). Check the database content.")
        else:
            QMessageBox.information(self, "palimpsestr",
                "SEF analysis complete. Loaded %d layer(s)." % n)

    def run_intrusions(self):
        if not self._check_provider():
            return
        path = self._require_sqlite()
        if not path:
            return
        import tempfile
        out = tempfile.mkdtemp(prefix="palimpsestr_")
        ip = os.path.join(out, "sef_intrusions.gpkg")
        params = {
            'Database_file': path,
            'Site': self._site(),
            'K': self.spin_k.value(),
            'Threshold': self.spin_thr.value(),
            'Intrusions': ip}
        try:
            res = processing.run(INTRUSIONS_ALG, params)
        except Exception as e:
            QMessageBox.critical(self, "palimpsestr", "Analysis failed:\n%s" % e)
            return
        n = self._load_outputs([(res.get('Intrusions', ip), "SEF intrusions", False)])
        if n == 0:
            QMessageBox.warning(self, "palimpsestr",
                "The analysis ran but produced no loadable layer.")
        else:
            QMessageBox.information(self, "palimpsestr",
                "Intrusion detection complete. Loaded %d layer(s)." % n)

    def _load_outputs(self, items):
        """Load (path, name, is_phase) outputs into the project; return count."""
        from qgis.core import QgsVectorLayer
        loaded = 0
        first = None
        for spec in items:
            pth, name, is_phase = spec
            if not pth or not os.path.exists(pth):
                continue
            lyr = QgsVectorLayer(pth, name, "ogr")
            if not lyr.isValid() or lyr.featureCount() == 0:
                continue
            QgsProject.instance().addMapLayer(lyr)
            loaded += 1
            if first is None:
                first = lyr
            if is_phase:
                self._style_phases(lyr)
        if first is not None:
            try:
                self.iface.setActiveLayer(first)
                self.iface.zoomToActiveLayer()
            except Exception:
                pass
        return loaded

    # --------------------------------------------------- install R scripts ---
    def _scripts_folder(self):
        try:
            from processing.tools.system import userFolder
            base = userFolder()
        except Exception:
            base = os.path.join(QgsApplication.qgisSettingsDirPath(), 'processing')
        folder = os.path.join(base, 'rscripts')
        os.makedirs(folder, exist_ok=True)
        return folder

    def install_scripts(self, silent=False):
        """Copy/overwrite the bundled .rsx into the Processing R scripts folder."""
        folder = self._scripts_folder()
        try:
            for name, content in RSX_SCRIPTS.items():
                with open(os.path.join(folder, name), 'w', encoding='utf-8') as f:
                    f.write(content)
            try:
                QgsApplication.processingRegistry().providerById('r').refreshAlgorithms()
            except Exception:
                pass
        except Exception as e:
            if not silent:
                QMessageBox.critical(self, "palimpsestr", "Could not install R scripts:\n%s" % e)
            return False
        if not silent:
            QMessageBox.information(self, "palimpsestr",
                "Installed/updated %d R scripts in:\n%s" % (len(RSX_SCRIPTS), folder))
        return True

    # ----------------------------------------------------------- styling ----
    def _style_phases(self, layer_id):
        layer = None
        if layer_id:
            layer = QgsProject.instance().mapLayer(layer_id) if isinstance(layer_id, str) else layer_id
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
#   self.actionPalimpsest = QAction("palimpsestr - Analisi palinsesti",
#                                    self.iface.mainWindow())
#   self.actionPalimpsest.triggered.connect(self.runPalimpsest)
#   self.analysisToolButton.addActions([self.actionPalimpsest])
#
#   def runPalimpsest(self):
#       from .tabs.Palimpsest import pyarchinit_Palimpsest
#       dlg = pyarchinit_Palimpsest(self.iface)
#       dlg.show()
#       self.pluginGui = dlg
# ---------------------------------------------------------------------------
