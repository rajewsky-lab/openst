import base64
import io

import numpy as np
from skimage.transform import SimilarityTransform

from openst.metadata.classes.base import BaseMetadata
from openst.utils.pimage import is_grayscale


class PairwiseAlignmentMetadata(BaseMetadata):
    def __init__(self, args):
        super().__init__(args)

        # it must be the same as the name of the HTML template under the templates dir
        self.metadata_type = "pairwise_alignment"
        self.alignment_results = []

    def add_alignment_result(self, alignment_result):
        self.alignment_results.append(alignment_result)

    def render(self):
        for alignment_result in self.alignment_results:
            alignment_result.render()


class AlignmentResult:
    def __init__(self, name, im_0, im_1, transformation_matrix, ransac_results, sift_results, keypoints0, keypoints1):
        self.name = name
        self.im_0 = im_0
        self.im_1 = im_1
        self.transformation_matrix = transformation_matrix
        self.ransac_results = ransac_results
        self.sift_results = sift_results
        self.keypoints0 = keypoints0
        self.keypoints1 = keypoints1
        self.alignment_rendered = None
        self.keypoints_rendered = None

    def visualize_alignment(self, fig=None, axes=None, show=False):
        import matplotlib.pyplot as plt
        from skimage.transform import warp

        # Apply the transformation matrix to the image
        if self.transformation_matrix is not None:
            aligned_im_1 = warp(self.im_1, SimilarityTransform(np.array(self.transformation_matrix)).inverse)
        else:
            aligned_im_1 = self.im_1

        # Display images before and after alignment
        if axes is None:
            fig, axes = plt.subplots(1, 2, figsize=(10, 5))
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

    def visualize_keypoints(self, fig=None, axes=None, show=False):
        import matplotlib.pyplot as plt
        from skimage.feature import plot_matches

        if axes is None:
            fig, axes = plt.subplots(1, 1)

        # TODO: automatically manage the number of channels to avoid errors here!
        if not is_grayscale(self.im_0):
            self.im_0 = self.im_0[..., 0]

        plot_matches(
            axes,
            self.im_0,
            self.im_1,
            self.keypoints0,
            self.keypoints1,
            np.repeat(np.arange(len(self.keypoints0)), 2).reshape(len(self.keypoints0), 2),
        ),

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
