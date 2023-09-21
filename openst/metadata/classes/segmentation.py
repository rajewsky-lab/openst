import base64
import io

import numpy as np

from openst.metadata.classes.base import BaseMetadata


class SegmentationMetadata(BaseMetadata):
    def __init__(self, args):
        super().__init__(args)

        # it must be the same as the name of the HTML template under the templates dir
        self.metadata_type = "pairwise_alignment"
        self.segmentation_result = None
        self.image_to_segment = None

    def add_alignment_result(self, alignment_result):
        self.alignment_results.append(alignment_result)

    def visualize_segmentation(self, fig=None, axes=None, show=False):
        import matplotlib.pyplot as plt

        # Display images before and after alignment
        if axes is None:
            fig, axes = plt.subplots(1, 2)
        axes[0].imshow(self.im_0, cmap="gray")
        axes[0].imshow(self.im_1, cmap="Blues_r", alpha=0.5)
        axes[0].set_axis_off()
        axes[0].set_title("Before alignment")
        axes[1].imshow(self.im_0, cmap="gray")
        axes[1].imshow(aligned_im_1, cmap="Blues_r", alpha=0.5)
        axes[1].set_axis_off()
        axes[1].set_title("After alignment")

        if show:
            plt.show()
        else:
            return fig, axes

    def plot_to_base64(self, fig):
        my_stringIObytes = io.BytesIO()
        fig.savefig(my_stringIObytes, format="jpg")
        my_stringIObytes.seek(0)
        return base64.b64encode(my_stringIObytes.read()).decode()

    def render(self):
        self.alignment_rendered, _ = self.visualize_alignment(axes=None, show=False)
        self.keypoints_rendered, _ = self.visualize_keypoints(axes=None, show=False)

        self.alignment_rendered = self.plot_to_base64(self.alignment_rendered)
        self.keypoints_rendered = self.plot_to_base64(self.keypoints_rendered)
