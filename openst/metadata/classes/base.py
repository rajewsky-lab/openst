import copy
import json
import logging
import sys

import numpy as np

from openst.utils.file import check_directory_exists

class BaseMetadata:
    def __init__(self, args):
        self.metadata_type = None
        self.cmdline = " ".join(sys.argv)
        self.args = args

    def render(self):
        """
        Render the content (not implemented).

        Raises:
            NotImplementedError: This function is not implemented.
        """
        raise NotImplementedError("This function is not implemented")

    def _get_dict_recursive(self, data):
        """
        Recursively convert an object and its attributes to a dictionary.

        Args:
            data (object): Input data object to be converted.

        Returns:
            dict or list: A dictionary representation of the input data or a
                          list of dictionaries for nested objects.
        """
        output_data = data
        if hasattr(data, "__dict__"):
            if not isinstance(data, dict):
                output_data = data.__dict__
            for k, v in output_data.items():
                output_data[k] = self._get_dict_recursive(v)
        elif isinstance(data, list):
            output_data = []
            for m in data:
                output_data.append(self._get_dict_recursive(m))
        elif isinstance(data, np.ndarray):
            logging.debug(
                """'np.ndarray' objects are not represented into a dictionary.
                        Hint: for large images, write a BaseMetadata.render method.
                            You can store compressed images as base64 strings. 
                            Other large datasets must not be stored as plain text."""
            )
            output_data = ""
        return output_data

    def copy(self):
        """
        Create a deep copy of the current object.

        Returns:
            object: A deep copy of the current object.
        """
        return copy.deepcopy(self)

    def get_dict(self):
        """
        Convert the current object and its attributes to a dictionary.

        Returns:
            dict or list: A dictionary representation of the object and its attributes.
        """
        data = self.copy()
        return self._get_dict_recursive(data)

    def save_json(self, path):
        """
        Save the object's data as JSON to a specified file path.

        Args:
            path (str): The file path where the JSON data will be saved.

        Raises:
            FileNotFoundError: If the specified directory in 'path' does not exist.
        """
        if not check_directory_exists(path):
            raise FileNotFoundError(f"A directory for {path} does not exist")

        json_content = self.get_dict()
        json_object = json.dumps(json_content, indent=4)

        with open(path, "w") as outfile:
            outfile.write(json_object)
