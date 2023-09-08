from html import HTML

# TODO: put these two classes into their own file
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
        from skimage.transform import SimilarityTransform

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

# TODO: import from its own

def load_pickle(file_path):
    """
    Load an object from a pickle file.

    Parameters:
        file_path (str): The path to the pickle file.

    Returns:
        any: The loaded object.
    """
    with open(file_path, "rb") as f:
        obj = pickle.load(f)
    return obj
        

def generate_alignment_report(alignment_result):
    alignment_report = HTML()

    # Add title
    alignment_report.h2("Alignment Result")

    # Add alignment parameters as a table
    alignment_params_table = HTML.table(
        [
            ["Parameter", "Value"],
            ["Image 0", alignment_result.im_0],
            ["Image 1", alignment_result.im_1],
            # Add more parameters as needed
        ],
        header_row=["Parameter", "Value"]
    )
    alignment_report.add(alignment_params_table)

    # Add images using visualize_alignment method
    # alignment_report.add(alignment_result.visualize_alignment())

    return str(alignment_report)

def generate_html_report(metadata):
    html_report = HTML()

    # Add title
    html_report.h1("Pairwise Alignment Report")

    # Add program parameters as a table
    program_params_table = HTML.table(
        [
            ["Parameter", "Value"],
            ["--image-in", metadata.args.image_in],
            ["--h5-in", metadata.args.h5_in],
            # Add more parameters as needed
        ],
        header_row=["Parameter", "Value"]
    )
    html_report.add(program_params_table)

    # Add alignment results
    for alignment_result in metadata.alignment_results:
        alignment_report = generate_alignment_report(alignment_result)
        html_report.add(alignment_report)

    return str(html_report)

if __name__ == "__main__":
    # Create a sample PairwiseAlignmentMetadata object for demonstration
    metadata = load_pickle("/data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/auto_aligned_he/fc_sts_063_4_alignment_metadata_loftr.pkl") # pass the pkl file here
    # Generate HTML report
    html_content = generate_html_report(metadata)
    # Save HTML report to a file with the name from the pkl
    with open("/data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/auto_aligned_he/fc_sts_063_4_alignment_report.html", "w") as f:
        f.write(html_content)