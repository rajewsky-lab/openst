from skimage.transform import SimilarityTransform


class PairwiseAlignmentMetadata:
    def __init__(self, args):
        self.args = args
        self.alignment_results = []

    def add_alignment_result(self, alignment_result):
        self.alignment_results.append(alignment_result)


class AlignmentResult:
    def __init__(self, im_0, im_1, transformation_matrix, ransac_results, sift_results, keypoints0, keypoints1):
        self.im_0 = im_0
        self.im_1 = im_1
        self.transformation_matrix = transformation_matrix
        self.ransac_results = ransac_results
        self.sift_results = sift_results
        self.keypoints0 = keypoints0
        self.keypoints1 = keypoints1

    def visualize_alignment(self, show=False):
        import matplotlib.pyplot as plt
        from skimage.transform import warp

        # Apply the transformation matrix to the image
        aligned_im_0 = warp(self.im_0, SimilarityTransform(self.transformation_matrix).inverse)

        # Display images before and after alignment
        _, axes = plt.subplots(1, 2)
        axes[0].imshow(self.im_0, cmap="gray")
        axes[0].imshow(self.im_1, cmap="Blues_r", alpha=0.5)
        axes[0].set_axis_off()
        axes[0].set_title("Before alignment")
        axes[1].imshow(aligned_im_0, cmap="gray")
        axes[1].imshow(self.im_1_after, cmap="Blues_r", alpha=0.5)
        axes[1].set_axis_off()
        axes[1].set_title("After alignment")

        if show:
            plt.show()
        else:
            return axes

    def visualize_keypoints(self, show=False):
        import matplotlib.pyplot as plt

        _, axes = plt.subplots(1, 2)
        axes[0].imshow(self.im_0, cmap="gray")
        axes[0].scatter(self.keypoints0[:, 0], self.keypoints0[:, 0])
        axes[0].set_title("Image 0")
        axes[1].imshow(self.im_1, cmap="gray")
        axes[1].scatter(self.keypoints1[:, 0], self.keypoints1[:, 0])
        axes[1].set_title("Image 1")

        if show:
            plt.show()
        else:
            return axes
