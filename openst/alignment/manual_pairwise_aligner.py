import json
import sys
import h5py
import gc
import signal
import numpy as np
import pyqtgraph as pg
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QHBoxLayout,
    QVBoxLayout,
    QFormLayout,
    QWidget,
    QPushButton,
    QFileDialog,
    QGraphicsEllipseItem,
    QSlider,
    QMessageBox,
    QGridLayout,
    QListView,
    QStyledItemDelegate,
    QLabel,
    QGroupBox,
    QProgressBar,
    QDialog,
    QCheckBox,
    QLineEdit,
    QScrollArea,
    QTreeView,
    QDialogButtonBox,
    QToolButton,
    QSizePolicy,
    QFrame,
)
from PyQt5.QtCore import (
    Qt,
    QThread,
    pyqtSignal,
    QParallelAnimationGroup,
    pyqtSlot,
    QAbstractAnimation,
    QPropertyAnimation,
)

from PyQt5.QtGui import QBrush, QColor, QStandardItemModel, QStandardItem, QIntValidator
from skimage.transform import warp

from openst.alignment.apply_transform import estimate_transform, apply_transform_to_coords, keypoints_json_to_dict
from openst.utils.pseudoimage import create_paired_pseudoimage
from openst.utils.file import h5_to_dict

# GUI elements
class CollapsibleBox(QWidget):
    def __init__(self, title="", parent=None):
        super(CollapsibleBox, self).__init__(parent)

        self.toggle_button = QToolButton(text=title, checkable=True, checked=False)
        self.toggle_button.setStyleSheet("QToolButton { border: none; }")
        self.toggle_button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.toggle_button.setArrowType(Qt.RightArrow)
        self.toggle_button.pressed.connect(self.on_pressed)

        self.toggle_animation = QParallelAnimationGroup(self)

        self.content_area = QScrollArea(maximumHeight=0, minimumHeight=0)
        self.content_area.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.content_area.setFrameShape(QFrame.NoFrame)

        lay = QVBoxLayout(self)
        lay.setSpacing(0)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.toggle_button)
        lay.addWidget(self.content_area)

        self.toggle_animation.addAnimation(QPropertyAnimation(self, b"minimumHeight"))
        self.toggle_animation.addAnimation(QPropertyAnimation(self, b"maximumHeight"))
        self.toggle_animation.addAnimation(QPropertyAnimation(self.content_area, b"maximumHeight"))

    @pyqtSlot()
    def on_pressed(self):
        checked = self.toggle_button.isChecked()
        self.toggle_button.setArrowType(Qt.DownArrow if not checked else Qt.RightArrow)
        self.toggle_animation.setDirection(QAbstractAnimation.Forward if not checked else QAbstractAnimation.Backward)
        self.toggle_animation.start()

    def setContentLayout(self, layout):
        lay = self.content_area.layout()
        del lay
        self.content_area.setLayout(layout)
        collapsed_height = self.sizeHint().height() - self.content_area.maximumHeight()
        content_height = layout.sizeHint().height()
        for i in range(self.toggle_animation.animationCount()):
            animation = self.toggle_animation.animationAt(i)
            animation.setDuration(500)
            animation.setStartValue(collapsed_height)
            animation.setEndValue(collapsed_height + content_height)

        content_animation = self.toggle_animation.animationAt(self.toggle_animation.animationCount() - 1)
        content_animation.setDuration(500)
        content_animation.setStartValue(0)
        content_animation.setEndValue(content_height)

class OverlayDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Please Wait")
        self.setWindowFlag(Qt.FramelessWindowHint)
        self.setAttribute(Qt.WA_TranslucentBackground)
        self.setStyleSheet("background-color: rgba(128, 128, 128, 128);")

        # Add a label with "Please wait" text
        self.label = QLabel("Please wait", self)
        self.label.setAlignment(Qt.AlignCenter)
        self.label.setStyleSheet("color: white; font-size: 20px;")
        self.progress_bar = QProgressBar(
            self,
            minimum=0,
            maximum=0,
            textVisible=False,
        )
        self.progress_bar.setAlignment(Qt.AlignCenter)

        # Create a layout for the label
        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.progress_bar)
        self.setLayout(layout)

    def updateTextLabel(self, text):
        self.label.setText(f"Please wait: {text}")

class TreeViewDialog(QDialog):
    def __init__(self, data, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select a Value")

        self.tree_view = QTreeView()
        self.tree_view.setModel(self.create_model(data))

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addWidget(self.tree_view)
        layout.addWidget(button_box)
        self.setLayout(layout)

    def create_model(self, data):
        # Implement a custom model to display the dictionary structure.
        # You can customize this part to match your dictionary format.
        model = QStandardItemModel()
        self.build_tree(data, model.invisibleRootItem())
        return model

    def build_tree(self, data, parent_item):
        if isinstance(data, dict):
            for key, value in data.items():
                item = QStandardItem(key)
                parent_item.appendRow(item)
                self.build_tree(value, item)
        elif isinstance(data, list):
            for value in data:
                self.build_tree(value, parent_item)

    def get_selected_path(self):
        # Get the selected path as a string.
        index = self.tree_view.selectionModel().currentIndex()
        if index.isValid():
            item = self.tree_view.model().itemFromIndex(index)
            path = []
            while item:
                path.insert(0, item.text())
                item = item.parent()
            return "/".join(path)

class InputDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.layername = QLineEdit(self)
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)

        layout = QFormLayout(self)
        layout.addRow("Layer name", self.layername)
        layout.addWidget(buttonBox)

        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)

    def getInputs(self):
        return self.layername.text()

class ItemDelegate(QStyledItemDelegate):
    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        option.decorationSize = option.rect.size()  # Adjust the size of items
        option.textElideMode = Qt.ElideNone  # Disable text eliding

class LabelledIntField(QWidget):
    def __init__(self, title, initial_value=None):
        QWidget.__init__(self)
        layout = QVBoxLayout()
        self.setLayout(layout)

        self.label = QLabel()
        self.label.setText(title)
        layout.addWidget(self.label)

        self.lineEdit = QLineEdit(self)
        self.lineEdit.setValidator(QIntValidator())
        if initial_value != None:
            self.lineEdit.setText(str(initial_value))
        layout.addWidget(self.lineEdit)
        layout.addStretch()

    def setLabelWidth(self, width):
        self.label.setFixedWidth(width)

    def setInputWidth(self, width):
        self.lineEdit.setFixedWidth(width)

    def getValue(self):
        return int(self.lineEdit.text())

# encoding numpy into json file
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

# image & pseudoimage renderer utils
class ImageRenderer(QThread):
    update_text = pyqtSignal(str)
    result_ready = pyqtSignal(dict)
    exception = pyqtSignal(Exception)

    def __init__(
        self,
        adata,
        recenter_coarse=False,
        rescale_factor_coarse=20,
        rescale_factor_fine=1,
        threshold_counts=1,
        pseudoimg_size=4000,
    ):
        super().__init__()
        self.adata = adata
        self.recenter_coarse = recenter_coarse
        self.rescale_factor_coarse = rescale_factor_coarse
        self.rescale_factor_fine = rescale_factor_fine
        self.threshold_counts = threshold_counts
        self.pseudoimg_size = pseudoimg_size
        self.spatial_path = "obs/spatial"
        self.img_path = "uns/spatial/image"
        self.layer = "all_tiles_coarse"

    def render_image_pair(
        self,
        in_coords: np.ndarray,
        total_counts: np.ndarray,
        tile_id: np.ndarray,
        staining_image: np.ndarray,
    ) -> dict:
        """
        Perform manual registration of spatial transcriptomics (STS) data with a staining image.

        Args:
            in_coords (np.ndarray): Input STS coordinates.
            total_counts (np.ndarray): Total UMI counts for each STS coordinate.
            tile_id (np.ndarray): Identifier for each STS coordinate. During the fine registration,
                    this 'tile_id' is used to aggregate the coordinates into buckets that
                    are aligned separately. Recommended for flow-cell based STS.
            staining_image (np.ndarray): Staining image for registration.

        Returns:
            image_pair (dict): Dictionary containing the pseudoimages, cropping and scaling limits, rendering scale.
        """
        image_pair = {}

        self.update_text.emit(f"Rendering '{self.layer}'")
        sts_coords = in_coords[:]
        _t_valid_coords = slice(None)

        if not self.recenter_coarse:
            _i_sts_coords_coarse_within_image_bounds = np.where(
                (sts_coords[:, 0] >= 0)
                & (sts_coords[:, 0] <= staining_image.shape[1])
                & (sts_coords[:, 1] >= 0)
                & (sts_coords[:, 1] <= staining_image.shape[0])
            )

            if len(_i_sts_coords_coarse_within_image_bounds[0]) < 10:
                raise ValueError("There are no spatial coordinates falling inside the image data.\n"+
                                 "Maybe spatial coordinates were never aligned to the image data.\n"+
                                 "Try enabling 'Is unaligned data' under 'Rendering settings' and render again.")
            sts_coords = sts_coords[_i_sts_coords_coarse_within_image_bounds]
            tile_id = tile_id[_i_sts_coords_coarse_within_image_bounds]
            total_counts = total_counts[_i_sts_coords_coarse_within_image_bounds]
            _t_valid_coords = tile_id == int(self.layer) if self.layer != 'all_tiles_coarse' else slice(None)

        min_lim, max_lim = sts_coords[_t_valid_coords].min(axis=0).astype(int), sts_coords[_t_valid_coords].max(axis=0).astype(int)
        x_all_min, y_all_min = min_lim
        x_all_max, y_all_max = max_lim

        # Preparing images and preprocessing routines
        _i_counts_above_threshold = total_counts > self.threshold_counts
        sts_coords = sts_coords[_i_counts_above_threshold]
        tile_id = tile_id[_i_counts_above_threshold]

        # We keep the limits also when filtering by UMI (avoid issues with offsets)
        sts_coords = np.concatenate([sts_coords, np.array([[x_all_max, y_all_max],
                                                        [x_all_min, y_all_min]])])
        
        if not self.recenter_coarse:
            min_lim, max_lim = [0, 0], staining_image.shape[:2][::-1]
            x_all_min, y_all_min = min_lim
            x_all_max, y_all_max = max_lim
            sts_coords = np.concatenate([sts_coords, np.array([[x_all_max, y_all_max],
                                                        [x_all_min, y_all_min]])])
    
        tile_id = np.concatenate([tile_id, [-1, -1, -2, -2]])

        # TODO: instead of this, we plot specific sections...
        if self.layer == "all_tiles_coarse":
            staining_image_rescaled = staining_image[:: self.rescale_factor_coarse, :: self.rescale_factor_coarse]
            sts_pseudoimage = create_paired_pseudoimage(
                sts_coords[:, ::-1],
                self.pseudoimg_size,
                staining_image_rescaled.shape,
                recenter=self.recenter_coarse,
                resize_method="cv2",
            )

            image_pair = {
                "imageA": staining_image_rescaled,
                "imageB": sts_pseudoimage["pseudoimage"],
                "lims": [0, 0, 0, 0], #[x_min, x_max, y_min, y_max],
                "factor_rescale": self.rescale_factor_coarse,
                "rescale_factor": sts_pseudoimage["rescale_factor"],
                "rescaling_factor": sts_pseudoimage["rescaling_factor"],
                "scale": self.pseudoimg_size,
                "offset_factor": sts_pseudoimage["offset_factor"]
            }
            self.update_text.emit(f"Creating metadata for '{self.layer}'")
        else:
            # Apply scaling to input image again, for fine registration
            staining_image_rescaled = staining_image[:: self.rescale_factor_fine, :: self.rescale_factor_fine]

            _t_valid_coords = (tile_id == int(self.layer)) | (tile_id == -1)

            _t_sts_pseudoimage = create_paired_pseudoimage(
                sts_coords[:, ::-1],  # we need to flip these coordinates
                self.pseudoimg_size,
                staining_image_rescaled.shape,
                _t_valid_coords,
                recenter=False,
                values=None,
                resize_method="cv2",
            )

            _t_sts_coords_to_transform = _t_sts_pseudoimage["coords_rescaled"] * _t_sts_pseudoimage["rescaling_factor"]

            # TODO: check if these can be from the original coordinates,
            # not from the offset (so the points are properly specified)
            min_lim, max_lim = _t_sts_coords_to_transform[_t_valid_coords].min(axis=0).astype(
                int
            ), _t_sts_coords_to_transform[_t_valid_coords].max(axis=0).astype(int)
            x_min, y_min = min_lim
            x_max, y_max = max_lim

            _pseudoimage = _t_sts_pseudoimage["pseudoimage"][x_min:x_max, y_min:y_max]

            self.update_text.emit(f"Creating metadata for '{self.layer}'")
            image_pair = {
                "imageA": staining_image_rescaled[x_min:x_max, y_min:y_max],
                "imageB": _pseudoimage,
                "lims": [x_min, x_max, y_min, y_max],
                "factor_rescale": self.rescale_factor_fine,
                "rescale_factor": _t_sts_pseudoimage["rescale_factor"],
                "rescaling_factor": _t_sts_pseudoimage["rescaling_factor"],
                "scale": self.pseudoimg_size,
                "offset_factor": _t_sts_pseudoimage["offset_factor"]
            }
            del _t_sts_pseudoimage
            gc.collect()

        self.update_text.emit(f"Done rendering '{self.layer}'")

        return image_pair

    def run(self):
        try:
            image_pair = self.render_image_pair(
                # put coordinates into XY (for correct rendering)
                self.adata[self.spatial_path][:][..., ::-1],
                self.adata["obs/total_counts"][:],
                self.adata["obs/tile_id/codes"][:],
                self.adata[self.img_path],
            )

            self.result_ready.emit(image_pair)
        except Exception as e:
            self.exception.emit(e)

class ColorImageView(pg.ImageView):
    """
    Wrapper around the ImageView to create a color lookup
    table automatically as there seem to be issues with displaying
    color images through pg.ImageView.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.lut = None

    def updateImage(self, autoHistogramRange=True):
        super().updateImage(autoHistogramRange)
        self.getImageItem().setLookupTable(self.lut)

# file management
class OpenAnndataWorkerThread(QThread):
    update_text = pyqtSignal(str)
    result_ready = pyqtSignal(h5py.File)

    def __init__(self, file_path):
        super().__init__()
        self.file_path = file_path

    def run(self):
        self.update_text.emit(f"Loading data")
        adata = self.load_image_pairs_from_h5ad()
        self.result_ready.emit(adata)

    def load_image_pairs_from_h5ad(self):
        return h5py.File(self.file_path, 'r+')
    
class OpenImageWorkerThread(QThread):
    update_text = pyqtSignal(str)
    result_ready = pyqtSignal(np.ndarray)

    def __init__(self, file_path):
        super().__init__()
        self.file_path = file_path

    def run(self):
        self.update_text.emit(f"Loading data")
        image = self.load_image_from_file()
        self.result_ready.emit(image)

    def load_image_from_file(self):
        # TODO: do lazy loading, so we don't need to load the whole image (also, support ome tiff etc)
        from tifffile import imread
        return imread(self.file_path)


class SavePointsWorkerThread(QThread):
    update_text = pyqtSignal(str)
    result_ready = pyqtSignal(bool)

    def __init__(self, file_path, points_to_write: dict):
        super().__init__()
        self.file_path = file_path
        self.points_to_write = points_to_write

    def run(self):
        self.update_text.emit(f"Saving points to text file")
        result = self.save_points_to_text_file()
        self.result_ready.emit(result)

    def save_points_to_text_file(self):
        # TODO: check if has extension, otherwise add   
        with open(f"{self.file_path}.json", "w") as fp:
            json.dump(self.points_to_write, fp, indent=4, cls=NpEncoder)

        return True


class ImageAlignmentApp(QMainWindow):
    def __init__(self, args):
        super().__init__()

        # Working variables
        self.imageA = None
        self.imageB = None
        self.mergedImage = None
        self.point_counter = 1
        self.image_pairs = {}
        self.active_view = None
        self.layer_names = []
        self.points_on_image = {}
        self.point_pairs = {}
        self.previous_layer = ""
        self.current_layer = ""
        self.points_to_write = ""
        self.adata = None
        self.adata_structure = None
        self.renderer = None
        self._merged_rgb_layer = None
        self._merged_pseudoimage_layer = None
        self.args = args

        # Initialize whole user interface
        self._init_ui()

        # Add all image pairs into layers
        for name in self.layer_names:
            self.add_image_pair(name)

        # Process arguments
        self._process_cmdline_args(args)

    # Building UI elements
    def _init_ui(self):
        self.worker_thread = None  # Keep track of the worker thread

        # UI elements
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)

        # Init UI elements
        self._init_layout()
        self._init_sidebar()
        self._init_imageviewers()

        # Set up keyPressEvent to handle backspace key
        self.keyPressEvent = self.on_key_press

    def _process_cmdline_args(self, args):
        if args.h5_in != "":
            self.open_h5ad(args.h5_in)
        
        if args.image_key != "":
            self._sidebar_buttons_tree_image.setText(args.image_key)

        if args.spatial_key != "":
            self._sidebar_buttons_tree_spatial.setText(args.spatial_key)

    def _init_layout(self):
        self.grid = QGridLayout(self.central_widget)

        self.grid.setColumnStretch(0, 1)
        self.grid.setColumnStretch(1, 1)

    def _init_sidebar(self):
        self._sidebar_hlayout = QHBoxLayout(self.central_widget)
        self._sidebar_layers_vlayout = QVBoxLayout(self.central_widget)
        self._sidebar_buttons_vlayout = QVBoxLayout(self.central_widget)
        self._sidebar_hlayout.addLayout(self._sidebar_layers_vlayout)
        self.grid.addLayout(self._sidebar_hlayout, 1, 1)

        # Add layers list
        self._init_layersviewer()

        # Controls group
        self._sidebar_buttons_groupbox = QGroupBox("Controls")
        self._sidebar_buttons_vlayout.addWidget(self._sidebar_buttons_groupbox)
        self._sidebar_buttons_groupbox_vbox = QVBoxLayout()
        self._sidebar_buttons_groupbox.setLayout(self._sidebar_buttons_groupbox_vbox)

        self._sidebar_scrollarea = QScrollArea()
        self._sidebar_scrollarea.setWidget(self._sidebar_buttons_groupbox)
        self._sidebar_scrollarea.setWidgetResizable(True)
        self._sidebar_hlayout.addWidget(self._sidebar_scrollarea)

        # Add a stretch
        self._sidebar_buttons_groupbox_vbox.addStretch()

        # Button: load STS
        self._sidebar_buttons_loadadata = QPushButton("Load h5ad", self)
        self._sidebar_buttons_loadadata.clicked.connect(self.open_h5ad)
        self._sidebar_buttons_groupbox_vbox.addWidget(self._sidebar_buttons_loadadata)

        # Properties of the h5ad file
        self._collapse_box_h5settings = CollapsibleBox("Data properties")
        lay = QVBoxLayout()
        self._sidebar_buttons_tree_image_label = QLabel("Image data", self)
        self._sidebar_buttons_tree_image = QPushButton("...", self)
        self._sidebar_buttons_tree_image.clicked.connect(self.select_image_path)
        lay.addWidget(self._sidebar_buttons_tree_image_label)
        lay.addWidget(self._sidebar_buttons_tree_image)

        self._sidebar_buttons_tree_load_image = QPushButton("Load image data", self)
        self._sidebar_buttons_tree_load_image.clicked.connect(self.open_image_data)
        lay.addWidget(self._sidebar_buttons_tree_load_image)

        self._sidebar_buttons_tree_spatial_label = QLabel("Spatial coordinates", self)
        self._sidebar_buttons_tree_spatial = QPushButton("...", self)
        self._sidebar_buttons_tree_spatial.clicked.connect(self.select_spatial_path)
        lay.addWidget(self._sidebar_buttons_tree_spatial_label)
        lay.addWidget(self._sidebar_buttons_tree_spatial)

        self._collapse_box_h5settings.setContentLayout(lay)
        self._sidebar_buttons_groupbox_vbox.addWidget(self._collapse_box_h5settings)

        # Fields: loading settings
        self._collapse_box_imagerender = CollapsibleBox("Rendering settings")
        lay = QVBoxLayout()
        self._sidebar_checkbox_offset_coarse = QCheckBox("Is unaligned data")
        self._sidebar_checkbox_offset_coarse.setChecked(False)
        self._sidebar_checkbox_offset_coarse.stateChanged.connect(self._update_imagerender_params)
        lay.addWidget(self._sidebar_checkbox_offset_coarse)

        self._sidebar_rescale_factor_fine = LabelledIntField("Rescale tiles", initial_value=1)
        self._sidebar_rescale_factor_fine.lineEdit.textChanged.connect(self._update_imagerender_params)
        lay.addWidget(self._sidebar_rescale_factor_fine)

        self._sidebar_rescale_factor_coarse = LabelledIntField("Rescale global", initial_value=20)
        self._sidebar_rescale_factor_coarse.lineEdit.textChanged.connect(self._update_imagerender_params)
        lay.addWidget(self._sidebar_rescale_factor_coarse)

        self._sidebar_threshold_counts = LabelledIntField("Threshold counts", initial_value=1)
        self._sidebar_threshold_counts.lineEdit.textChanged.connect(self._update_imagerender_params)
        lay.addWidget(self._sidebar_threshold_counts)

        self._sidebar_pseudoimage_size = LabelledIntField("Pseudoimage size", initial_value=4000)
        self._sidebar_pseudoimage_size.lineEdit.textChanged.connect(self._update_imagerender_params)
        lay.addWidget(self._sidebar_pseudoimage_size)

        self.update_opacity_A_slider_label = QLabel("Image A opacity (merge)", self)
        self.update_opacity_A_slider = QSlider(Qt.Horizontal)
        self.update_opacity_A_slider.setMinimum(0)
        self.update_opacity_A_slider.setMaximum(100)
        self.update_opacity_A_slider.setValue(100)
        self.update_opacity_A_slider.valueChanged.connect(self.update_opacity_A)
        lay.addWidget(self.update_opacity_A_slider_label)
        lay.addWidget(self.update_opacity_A_slider)

        self.update_opacity_B_slider_label = QLabel("Image B opacity (merge)", self)
        self.update_opacity_B_slider = QSlider(Qt.Horizontal)
        self.update_opacity_B_slider.setMinimum(0)
        self.update_opacity_B_slider.setMaximum(100)
        self.update_opacity_B_slider.setValue(50)
        self.update_opacity_B_slider.valueChanged.connect(self.update_opacity_B)
        lay.addWidget(self.update_opacity_B_slider_label)
        lay.addWidget(self.update_opacity_B_slider)

        self._collapse_box_imagerender.setContentLayout(lay)
        self._sidebar_buttons_groupbox_vbox.addWidget(self._collapse_box_imagerender)

        # Section: keypoints
        self._collapse_box_keypoints = CollapsibleBox("Keypoints properties")
        lay = QVBoxLayout()
        self._sidebar_buttons_loadbutton = QPushButton("Load keypoints", self)
        self._sidebar_buttons_loadbutton.clicked.connect(self.load_point_pairs)
        lay.addWidget(self._sidebar_buttons_loadbutton)

        self._sidebar_buttons_savebutton = QPushButton("Save keypoints", self)
        self._sidebar_buttons_savebutton.clicked.connect(self.save_point_pairs)
        lay.addWidget(self._sidebar_buttons_savebutton)

        self.point_size_slider_label = QLabel("Point size", self)
        self.point_size_slider = QSlider(Qt.Horizontal)
        self.point_size_slider.setMinimum(1)
        self.point_size_slider.setMaximum(200)
        self.point_size_slider.setValue(100)
        lay.addWidget(self.point_size_slider_label)
        lay.addWidget(self.point_size_slider)

        self._collapse_box_keypoints.setContentLayout(lay)
        self._sidebar_buttons_groupbox_vbox.addWidget(self._collapse_box_keypoints)

        # Button: rendering
        self._sidebar_buttons_render = QPushButton("Render", self)
        self._sidebar_buttons_render.clicked.connect(self.render)
        self._sidebar_buttons_groupbox_vbox.addWidget(self._sidebar_buttons_render)

        # Button: align merge
        self._sidebar_buttons_preview_alignment = QPushButton("Preview alignment", self)
        self._sidebar_buttons_preview_alignment.clicked.connect(self.preview_alignment)
        self._sidebar_buttons_groupbox_vbox.addWidget(self._sidebar_buttons_preview_alignment)

        # Button: apply transform to data
        self._sidebar_buttons_apply_to_data = QPushButton("Apply to data", self)
        self._sidebar_buttons_apply_to_data.clicked.connect(self.apply_to_data)
        self._sidebar_buttons_groupbox_vbox.addWidget(self._sidebar_buttons_apply_to_data)

    def select_image_path(self):
        self.show_path_selection_dialog(self._sidebar_buttons_tree_image)

    def select_spatial_path(self):
        self.show_path_selection_dialog(self._sidebar_buttons_tree_spatial)

    def show_path_selection_dialog(self, button):
        if self.adata is None:
            QMessageBox.warning(
                self,
                "Warning",
                "First, load a valid AnnData h5 file!",
            )
            return False

        dialog = TreeViewDialog(self.adata_structure, self)
        if dialog.exec_() == QDialog.Accepted:
            selected_path = dialog.get_selected_path()
            button.setText(selected_path)
            self._update_imagerender_params()
            # and set the path to the h5file

        return True

    def _init_layersviewer(self):
        self._sidebar_layers_groupbox = QGroupBox("Layer selector")
        self._sidebar_layers_vlayout.addWidget(self._sidebar_layers_groupbox)

        self._sidebar_layers_groupbox_vbox = QVBoxLayout()
        self._sidebar_layers_groupbox.setLayout(self._sidebar_layers_groupbox_vbox)
        self._sidebar_layers_label = QLabel("Layers")
        self._sidebar_layers_listview = QListView(self)
        self._sidebar_layers_listview.setEditTriggers(QListView.NoEditTriggers)
        self._sidebar_layers_groupbox_vbox.addWidget(self._sidebar_layers_listview)
        delegate = ItemDelegate()
        self._sidebar_layers_listview.setItemDelegate(delegate)

        # Data model
        self._sidebar_layers_listmodel = QStandardItemModel()

        # Connecting models and UI changes
        self._sidebar_layers_listview.setModel(self._sidebar_layers_listmodel)
        self._sidebar_layers_listview.selectionModel().selectionChanged.connect(self.update_selected_layer)

    def _init_imageviewers(self):
        def _hide_ui_elements(imageview):
            # imageview.ui.histogram.hide()
            imageview.ui.roiBtn.hide()
            imageview.ui.menuBtn.hide()

        # Add viewers for the images A and B (source, dest)
        self.image_view_a = ColorImageView(self)
        self.grid.addWidget(self.image_view_a, 0, 0)
        _hide_ui_elements(self.image_view_a)

        self.image_view_b = ColorImageView(self)
        self.grid.addWidget(self.image_view_b, 0, 1)
        _hide_ui_elements(self.image_view_b)

        self.image_view_a.mouseDoubleClickEvent = lambda e: self.set_active_view(e, "A")
        self.image_view_b.mouseDoubleClickEvent = lambda e: self.set_active_view(e, "B")

        # Add viewers for the merged image after alignment
        self.image_view_merged = ColorImageView(self)
        self.grid.addWidget(self.image_view_merged, 1, 0)
        _hide_ui_elements(self.image_view_merged)

        # Synchronize axes of the image viewers
        self.syncedPlots = [
            self.image_view_a,
            self.image_view_b,
            self.image_view_merged,
        ]
        self.syncedPlots[0].view.setXLink(self.syncedPlots[1].view)
        self.syncedPlots[1].view.setXLink(self.syncedPlots[2].view)
        self.syncedPlots[0].view.setYLink(self.syncedPlots[1].view)
        self.syncedPlots[1].view.setYLink(self.syncedPlots[2].view)

    def preview_alignment(self):
        # get keypoints and estimate transformation matrix
        _current_keypoints = self.points_on_image[self.current_layer]
        _current_point_pairs = self.point_pairs[self.current_layer]

        if len(_current_keypoints) < 4 or len(_current_keypoints) % 2 != 0:
            QMessageBox.warning(self, "Warning", "Needs >= 2 keypoints to estimate alignment model")
            return

        _t_mkpts0 = np.zeros((int(len(_current_keypoints)/2), 2))
        _t_mkpts1 = np.zeros((int(len(_current_keypoints)/2), 2))

        for i, keypoint_i in enumerate(range(0, len(_current_keypoints), 2)):
            pointA, pointB = _current_keypoints[keypoint_i], _current_keypoints[keypoint_i + 1]
    
            xA_0, yA_0 = _current_point_pairs[keypoint_i]
            xB_0, yB_0 = _current_point_pairs[keypoint_i + 1]

            bounding_rect_A = pointA.boundingRect()
            bounding_rect_B = pointB.boundingRect()

            width_A = bounding_rect_A.width()
            radius_A = width_A/2
            width_B = bounding_rect_B.width()
            radius_B = width_B/2

            xA, yA = pointA.pos().x(), pointA.pos().y()
            xB, yB = pointB.pos().x(), pointB.pos().y()
            _t_mkpts0[i] = np.array([xA+xA_0+radius_A, yA+yA_0+radius_A])
            _t_mkpts1[i] = np.array([xB+xB_0+radius_B, yB+yB_0+radius_B])

        _t_image_B = self.imageB
        _t_matrix, needs_flip = estimate_transform("similarity", _t_mkpts1[:, ::-1], _t_mkpts0[:, ::-1])
        if needs_flip:
            _t_image_B = _t_image_B[::-1]

        # warp the image using the transformation matrix (similarity)
        _t_image_B = warp(_t_image_B, _t_matrix.inverse) # warp from skimage

        # Remove current images from merged viewport (if any)
        if self._merged_rgb_layer is not None:
            self.image_view_merged.removeItem(self._merged_rgb_layer)
            self.image_view_merged.removeItem(self._merged_pseudoimage_layer)
            self.image_view_merged.clear()
        
        # Render merged
        self._merged_rgb_layer = pg.ImageItem(self.imageA)
        self._merged_rgb_layer.setOpts(update=True, opacity=self.update_opacity_A_slider.value()/100)
        self.image_view_merged.addItem(self._merged_rgb_layer)
        self._merged_pseudoimage_layer = pg.ImageItem(_t_image_B)
        self._merged_pseudoimage_layer.setOpts(update=True, opacity=self.update_opacity_B_slider.value()/100)
        self.image_view_merged.addItem(self._merged_pseudoimage_layer)

    def apply_to_data(self):
        # so we can validate that applying the transform will look good...
        self.preview_alignment()

        # we show a dialog to get the name for the new layer
        # TOOD: get dialog
        dialog = InputDialog()
        if dialog.exec():
            spatial_key_out = dialog.getInputs()
        else:
            QMessageBox.warning(self, "Warning", "Please specify a name to save the coordinates into 'obsm'\n(e.g., 'obsm/spatial_manual_coarse')")
            return
    
        # we apply the transformation using the manual_pairwise_aligner code
        in_coords = self.adata[self.renderer.spatial_path][:]

        # TODO: we only do for the whole display layer, we need to adapt for the whole thing...
        tile_id = None
        keypoints = keypoints_json_to_dict(self.generate_point_pairs_dict()['points'])

        # Apply transform to all coordinates & retransform back
        sts_coords_coarse = in_coords[..., ::-1].copy()
        sts_transformed = apply_transform_to_coords(sts_coords_coarse, tile_id, keypoints, check_bounds=False)

        if spatial_key_out in self.adata:
            self.adata[spatial_key_out][...] = sts_transformed[:][..., ::-1]
        else:
            self.adata[spatial_key_out] = sts_transformed[:][..., ::-1]

        # TODO: display some success feedback
        self.adata_structure = h5_to_dict(self.adata)

    def _update_imagerender_params(self):
        if self.adata is None or self.renderer is None:
            QMessageBox.warning(
                self,
                "Warning",
                "You must load a valid AnnData h5 file for these parameters to have an effect",
            )

        self.renderer.recenter_coarse = self._sidebar_checkbox_offset_coarse.isChecked()
        self.renderer.rescale_factor_coarse = int(self._sidebar_rescale_factor_coarse.getValue())
        self.renderer.rescale_factor_fine = int(self._sidebar_rescale_factor_fine.getValue())
        self.renderer.threshold_counts = int(self._sidebar_threshold_counts.getValue())
        self.renderer.pseudoimg_size = int(self._sidebar_pseudoimage_size.getValue())
        self.renderer.spatial_path = self._sidebar_buttons_tree_spatial.text()
        self.renderer.img_path = self._sidebar_buttons_tree_image.text()

    # Interactive data loading
    def prepare_layers(self, layers):
        self.layer_names = list(layers)
        for name in self.layer_names:
            self.add_image_pair(name)

        self.points_on_image = {i: [] for i in self.layer_names}
        self.point_pairs = {i: [] for i in self.layer_names}

        self.current_layer = self.layer_names[0]
        self.previous_layer = self.current_layer

    def data_loaded(self, adata):
        self.adata = adata

        # Get the layer names from the tile_id
        layers = ["all_tiles_coarse"]
        layers += [f"{i}" for i in np.unique(self.adata["obs/tile_id/codes"])]
        self.prepare_layers(layers)

        # Setup the tree structure
        self.adata_structure = h5_to_dict(self.adata)

        self.renderer = ImageRenderer(self.adata)
        self.renderer.update_text.connect(self.overlay_dialog.updateTextLabel)
        self.renderer.finished.connect(self.overlay_dialog.accept)
        self.renderer.result_ready.connect(self.display_images)
        self.renderer.exception.connect(self.handle_render_exception)
        self._update_imagerender_params()

    def image_data_loaded(self, img):
        # TODO: check if it is a dictionary or what, so we can save accordingly
        if self.adata is None:
            QMessageBox.warning(self, "Warning", "No anndata file was loaded.")
            return
    
        dialog = InputDialog()
        if dialog.exec():
            img_key_out = dialog.getInputs()
        else:
            QMessageBox.warning(self, "Warning", "Please specify a path to save the image into 'uns'")
            return

        # TODO: do the copy lazily into the h5 file (no need to load the whole image into memory)
        if f"uns/{img_key_out}" in self.adata:
            self.adata[f"uns/{img_key_out}"][...] = img
        else:
            self.adata[f"uns/{img_key_out}"] = img
    
        # Setup the tree structure
        self.adata_structure = h5_to_dict(self.adata)

    def open_h5ad(self, file_path=False):
        if file_path is False:
            options = QFileDialog.Options()
            file_path, _ = QFileDialog.getOpenFileName(self, "Open File", "", "AnnData h5 files (*.h5ad)", options=options)

        if file_path:
            self.overlay_dialog = OverlayDialog(self)
            self.overlay_dialog.setWindowModality(Qt.WindowModal)
            self.overlay_dialog.show()

            self.worker_thread = OpenAnndataWorkerThread(file_path)
            self.worker_thread.result_ready.connect(self.data_loaded)

            self.worker_thread.update_text.connect(self.overlay_dialog.updateTextLabel)
            self.worker_thread.finished.connect(self.overlay_dialog.accept)

            self.worker_thread.start()

    def open_image_data(self, file_path=False):
        if file_path is False:
            options = QFileDialog.Options()
            file_path, _ = QFileDialog.getOpenFileName(self, "Open File", "", "Image files (*.tiff *.tif *.png *.jpeg *.ome.tiff)", options=options)

        if file_path:
            self.overlay_dialog = OverlayDialog(self)
            self.overlay_dialog.setWindowModality(Qt.WindowModal)
            self.overlay_dialog.show()

            self.worker_thread = OpenImageWorkerThread(file_path)
            self.worker_thread.result_ready.connect(self.image_data_loaded)

            self.worker_thread.update_text.connect(self.overlay_dialog.updateTextLabel)
            self.worker_thread.finished.connect(self.overlay_dialog.accept)

            self.worker_thread.start()

    # Interactive data displaying
    def set_active_view(self, event, which):
        self.active_view = which
        self.mouseDoubleClickEvent(event)

    def update_selected_layer(self, selected, deselected):
        selected_item = self._sidebar_layers_listmodel.itemFromIndex(selected.indexes()[0])
        item_text = selected_item.text()

        self.previous_layer = self.current_layer
        self.current_layer = item_text
        self.renderer.layer = self.current_layer
        # self.renderer.run()

    def render(self):
        self.overlay_dialog = OverlayDialog(self)
        self.overlay_dialog.setWindowModality(Qt.WindowModal)
        self.overlay_dialog.show()

        self.renderer.update_text.connect(self.overlay_dialog.updateTextLabel)
        self.renderer.finished.connect(self.overlay_dialog.accept)

        self.renderer.start()

    def handle_render_exception(self, exception):
        QMessageBox.warning(self, "Rendering exception", str(exception))

    def display_images(self, image_pair=None):
        if image_pair is None:
            return

        self.image_pairs[self.current_layer] = image_pair

        self.imageA = image_pair["imageA"]
        self.imageB = image_pair["imageB"]

        # Show image A
        self.image_view_a.setImage(image_pair["imageA"])
        self.image_view_a.getView().invertY(True)
        self.image_view_a.getView().setAspectLocked(True)

        # Show image B
        self.image_view_b.setImage(image_pair["imageB"])
        self.image_view_b.getView().invertY(True)
        self.image_view_b.getView().setAspectLocked(True)

        # Remove current images from merged viewport (if any)
        if self._merged_rgb_layer is not None:
            self.image_view_merged.removeItem(self._merged_rgb_layer)
            self.image_view_merged.removeItem(self._merged_pseudoimage_layer)
            self.image_view_merged.clear()
    
        # Show image merged
        self._merged_rgb_layer = pg.ImageItem(self.imageA)
        self._merged_rgb_layer.setOpts(update=True, opacity=self.update_opacity_A_slider.value()/100)
        self.image_view_merged.addItem(self._merged_rgb_layer)
        self._merged_pseudoimage_layer = pg.ImageItem(self.imageB)
        self._merged_pseudoimage_layer.setOpts(update=True, opacity=self.update_opacity_B_slider.value()/100)
        self.image_view_merged.addItem(self._merged_pseudoimage_layer)

        # self.image_view_merged.setImage(self.mergedImage)
        # self.image_view_merged.getView().invertY(True)
        # self.image_view_merged.getView().setAspectLocked(True)

        # Clear existing points on the image
        if self.previous_layer != "":
            for point in self.points_on_image[self.previous_layer]:
                self.image_view_a.getView().removeItem(point)
                self.image_view_b.getView().removeItem(point)

        # Remove the points added to the ImageView object
        for _, v in self.points_on_image.items():
            for o in v:
                try:
                    self.image_view_a.removeItem(o)
                except:
                    pass
                
                try:
                    self.image_view_b.removeItem(o)
                except:
                    pass
        
        # Add points previously stored in the layer
        if self.current_layer != "":
            _points = self.points_on_image[self.current_layer].copy()
            del self.points_on_image[self.current_layer]
            self.points_on_image[self.current_layer] = []
            for i, _point in enumerate(_points):
                if i % 2 == 0:
                    self.image_view_a.getView().addItem(_point)
                else:
                    self.image_view_b.getView().addItem(_point)
                self.points_on_image[self.current_layer].append(_point)

        # Reset images from the object, since they will not be used for anything else
        self.image_pairs[self.current_layer]["imageA"] = None
        self.image_pairs[self.current_layer]["imageB"] = None

    def update_opacity_A(self, opacity_value):
        # The opacity value needs to be between 0-100 because QSlider does not support float
        self._merged_rgb_layer.setOpts(opacity=opacity_value/100)
        self.image_view_merged.update()

    def update_opacity_B(self, opacity_value):
        self._merged_pseudoimage_layer.setOpts(opacity=opacity_value/100)
        self.image_view_merged.update()

    def load_point_pairs_dict(self, point_pairs_dict):
        for keypoint in point_pairs_dict['points']:
            # set the layer and add the point
            self.current_layer = keypoint['layer']

            # readjust with the point size (will be added later)
            _p_sz = (self.point_size_slider.value() / 2)
            self.add_point((float(keypoint['point_src'][0]) + _p_sz, float(keypoint['point_src'][1]) + _p_sz), which="A")
            self.add_point((float(keypoint['point_dst'][0]) + _p_sz, float(keypoint['point_dst'][1]) + _p_sz), which="B")

    def load_point_pairs(self):
        options = QFileDialog.Options()
        filedialog = QFileDialog()
        filedialog.setDefaultSuffix(".json")
        filedialog.setNameFilters(["JSON (*.json)"])
        file_path, _ = filedialog.getOpenFileName(self, "Load File", "", "JSON (*.json)", options=options)

        def load_keypoints_from_json(fname: str):
            with open(fname) as j:
                _keypoints_dict = json.load(j)

            return _keypoints_dict
        
        points_to_load = load_keypoints_from_json(file_path)
        _old_layer = self.current_layer
        self.load_point_pairs_dict(points_to_load)
        self.current_layer = _old_layer

    def generate_point_pairs_dict(self):
        points_to_write = {"points": []}
        for k, v in self.image_pairs.items():
            for i in range(0, len(self.points_on_image[k]), 2):
                pointA = self.points_on_image[k][i]
                pointB = self.points_on_image[k][i + 1]

                bounding_rect_A = pointA.boundingRect()
                bounding_rect_B = pointB.boundingRect()

                width_A = bounding_rect_A.width()
                radius_A = width_A/2
                width_B = bounding_rect_B.width()
                radius_B = width_B/2

                xA_0, yA_0 = self.point_pairs[k][i]
                xB_0, yB_0 = self.point_pairs[k][i + 1]

                xA, yA = pointA.pos().x(), pointA.pos().y()
                xB, yB = pointB.pos().x(), pointB.pos().y()

                _ofs_factor = v["offset_factor"]
                _ofs_factor = np.array([0, 0]) if _ofs_factor is None else _ofs_factor
                _rescale = v["factor_rescale"]
                _lims = v["lims"]

                points_to_write["points"].append(
                    {
                        "layer": k,
                        "point_src": [f"{(xA+xA_0+radius_A):.2f}", f"{(yA+yA_0+radius_A):.2f}"],
                        "point_dst": [f"{(xB+xB_0+radius_B):.2f}", f"{(yB+yB_0+radius_B):.2f}"],
                        "lims": _lims,
                        "factor_rescale": _rescale,
                        "offset_factor": _ofs_factor,
                        "rescaling_factor": v['rescaling_factor'],
                        "scale": v['scale'],
                        "rescale_factor": v['rescale_factor'],
                        "point_src_offset_rescaled": [
                            f"{((((xA+xA_0+radius_A)+_lims[0])*_rescale)):.2f}",
                            f"{((((yA+yA_0+radius_A)+_lims[2])*_rescale)):.2f}",
                        ],
                        "point_dst_offset_rescaled": [
                            # f"{((((xB+xB_0)+_lims[0])*_rescale*_scale)+_ofs_factor[1]):.2f}",
                            # f"{((((yB+yB_0)+_lims[2])*_rescale*_scale)+_ofs_factor[0]):.2f}",
                            f"{(((xB+xB_0+radius_B)/(v['rescaling_factor']*v['scale']))*v['rescale_factor']+_lims[0]+_ofs_factor[0]):.2f}",
                            f"{(((yB+yB_0+radius_B)/(v['rescaling_factor']*v['scale']))*v['rescale_factor']+_lims[2]+_ofs_factor[1]):.2f}"
                        ],
                    }
                )

        return points_to_write
    
    def save_point_pairs(self):
        if all(len(value) == 0 for value in self.points_on_image.values()):
            QMessageBox.warning(self, "Warning", "No points to save.")
            return

        points_to_write = self.generate_point_pairs_dict()

        options = QFileDialog.Options()
        filedialog = QFileDialog()
        filedialog.setDefaultSuffix(".json")
        filedialog.setNameFilters(["JSON (*.json)"])
        file_path, _ = filedialog.getSaveFileName(self, "Save File", "", "JSON (*.json)", options=options)

        self.overlay_dialog = OverlayDialog(self)
        self.overlay_dialog.setWindowModality(Qt.WindowModal)
        self.overlay_dialog.show()
        self.worker_thread = SavePointsWorkerThread(file_path, points_to_write)
        self.worker_thread.finished.connect(self.overlay_dialog.accept)
        self.worker_thread.start()

    def add_image_pair(self, name):
        self._sidebar_layers_listmodel.appendRow(QStandardItem(name))

    def add_point_to_image(self, pos, which="B"):
        x, y = pos
        point_size = self.point_size_slider.value()
        point = QGraphicsEllipseItem(x - (point_size / 2), y - (point_size / 2), point_size, point_size)
        if which == "A":
            point.setBrush(QBrush(QColor(0, 0, 255, 128)))
            point.setFlag(QGraphicsEllipseItem.ItemIsMovable)
            self.image_view_a.getView().addItem(point)
        elif which == "B":
            point.setBrush(QBrush(QColor(255, 0, 0, 128)))
            point.setFlag(QGraphicsEllipseItem.ItemIsMovable)
            self.image_view_b.getView().addItem(point)
            self.point_counter += 1

        self.points_on_image[self.current_layer].append(point)

    def add_point(self, pos, which="B"):
        x, y = pos
        point_size = self.point_size_slider.value()
        self.add_point_to_image(pos, which)
        self.point_pairs[self.current_layer].append((x - (point_size / 2), y - (point_size / 2)))

    def on_key_press(self, event):
        if event.key() == Qt.Key_Backspace and len(self.points_on_image[self.current_layer]) > 0:
            last_point = self.points_on_image[self.current_layer].pop()
            if len(self.point_pairs[self.current_layer]) % 2 == 0:
                self.image_view_b.getView().removeItem(last_point)
            else:
                self.image_view_a.getView().removeItem(last_point)
            self.point_pairs[self.current_layer].pop()

    def mouseDoubleClickEvent(self, event):
        if self.imageA is None or self.imageB is None:
            return

        if self.active_view is None:
            return

        if self.active_view == "A" and len(self.point_pairs[self.current_layer]) % 2 == 0:
            point_img_coord = self.image_view_a.getImageItem().mapFromScene(
                self.image_view_a.getView().mapToScene(event.pos())
            )
            self.add_point((point_img_coord.x(), point_img_coord.y()), which="A")
        elif self.active_view == "B" and len(self.point_pairs[self.current_layer]) % 2 != 0:
            point_img_coord = self.image_view_b.getImageItem().mapFromScene(
                self.image_view_b.getView().mapToScene(event.pos())
            )
            self.add_point((point_img_coord.x(), point_img_coord.y()), which="B")
        elif self.active_view == "B":
            QMessageBox.warning(
                self,
                "Warning",
                "First add a new keypoint to A (left)",
            )
        elif self.active_view == "A":
            QMessageBox.warning(
                self,
                "Warning",
                "First add a corresponding keypoint to the image B (right)",
            )


def _run_manual_pairwise_aligner(args):
    app = QApplication(sys.argv)
    window = ImageAlignmentApp(args)
    window.setWindowTitle("Manual Pairwise Alignment (Open-ST)")
    window.setGeometry(100, 100, 1000, 1000)
    window.show()

    # Catch KeyboardInterrupt (Ctrl+C) and close the application gracefully
    def handle_interrupt(signal, frame):
        if window.worker_thread and window.worker_thread.isRunning():
            window.worker_thread.terminate()  # Terminate the worker thread if running
        app.quit()  # Quit the application gracefully

    signal.signal(signal.SIGINT, handle_interrupt)  # Register the signal handler

    sys.exit(app.exec_())

if __name__ == "__main__":
    from openst.cli import get_manual_pairwise_aligner_parser
    args = get_manual_pairwise_aligner_parser().parse_args()
    _run_manual_pairwise_aligner(args)