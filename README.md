[![Downloads](https://pepy.tech/badge/openst)](https://pepy.tech/project/openst)
[![PyPI Version](https://img.shields.io/pypi/v/openst.svg)](https://pypi.org/project/openst)
[![PyPI License](https://img.shields.io/pypi/l/openst.svg)](https://pypi.org/project/openst)

<img
  src="https://raw.githubusercontent.com/rajewsky-lab/openst/main/docs/static/img/openst_workflow.svg"
  class="dark-light" align="right" width="350" alt="image"
/>

# Open-ST: democratizing spatial transcriptomics
Open-ST is an open-source [spatial transcriptomics](https://en.wikipedia.org/wiki/Spatial_transcriptomics) method 
with efficient whole-transcriptome capture at sub-cellular resolution (0.5 µm) at low cost 
(<150 Euro library preparation per 12 mm^2).

Open-ST requires standard lab equipment
and includes open-source software for seamless data processing and analysis.

- Discuss development on [GitHub](https://github.com/rajewsky-lab/openst).
- Read the [step-by-step protocol](https://openst.github.io).
- Install via `pip install openst`.

## Tutorials and examples
We provide [various open-ST datasets](https://openst.github.io/examples) collected by the [Rajewsky lab @ MDC Berlin](https://www.mdc-berlin.de/n-rajewsky), from various tissues/organisms.
There, links to raw and processed data are available, as well as step-by-step guides for (pre)processing and downstream analysis. 

## Documentation
All the detail to the experimental and computational aspects of our method are available in [our website](https://openst.github.io).

We love to have an open approach to documentation. We decided to use [mkdocs](https://github.com/mkdocs/mkdocs) as our documentation backend 
to make your life easier. So, feel free to suggest changes by opening a 
[documentation-related issue](https://github.com/rajewsky-lab/openst/issues/new?assignees=&labels=docs&template=&title=)!

You can build the documentation locally by following these steps:
1. Clone this repository
  ```sh
  git clone https://github.com/rajewsky-lab/openst
  ```
2. Install the dependencies for building the documentation:
   ```sh
   pip install mkdocs-material mkdocs-autorefs
   ```
3. Serve `mkdocs` with the following command:
   ```sh
   mkdocs serve -f openst/mkdocs.yml
   ```

## Contributing
Open-ST is an open-source project mostly maintained by the [Rajewsky lab @ MDC Berlin](https://www.mdc-berlin.de/n-rajewsky) - so, your involvement is warmly welcome! 
If you're excited to join us, we recommend the following steps:

- Looking for ideas? See our [Volunteer Project Board](https://github.com/orgs/rajewsky-lab/projects/1) to see what we may need help with.
- Found a bug? Contact an admin in the form of an [issue](https://github.com/rajewsky-lab/openst/issues/new?assignees=&labels=&template=bug-report.md&title=).
- Implement your idea following guidelines set by the [official contributing guide](CONTRIBUTING.md)
- Wait for admin approval; approval is iterative, but if accepted will belong to the main repository.

In general, you can always refer to the [contribution guidelines](CONTRIBUTING.md) for more details!
Currently, only [admins](https://github.com/orgs/rajewsky-lab/people) will be merging all accepted changes.

## License
The software tools of this project are under the GNU License - see the [LICENSE](LICENSE) file for details.