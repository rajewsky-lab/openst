[![docs](https://github.com/rajewsky-lab/openst/actions/workflows/docs.yml/badge.svg)](https://github.com/rajewsky-lab/openst/actions/workflows/docs.yml)
[![Publish Docker image](https://github.com/rajewsky-lab/openst/actions/workflows/docker_hub.yml/badge.svg)](https://github.com/rajewsky-lab/openst/actions/workflows/docker_hub.yml)
[![Downloads](https://pepy.tech/badge/openst)](https://pepy.tech/project/openst)
[![PyPI Version](https://img.shields.io/pypi/v/openst.svg)](https://pypi.org/project/openst)
[![PyPI License](https://img.shields.io/pypi/l/openst.svg)](https://pypi.org/project/openst)
[![Gitter chat](https://badges.gitter.im/openst/community.png)](https://gitter.im/openst/community)

<img
  src="https://raw.githubusercontent.com/rajewsky-lab/openst/5617df9d35341778487d4c623eaac53b61000006/docs/static/img/openst_logo_color.png"
  class="dark-light" align="right" width="350" alt="image"
/>

# Open-ST: open-source spatial transcriptomics

### [🌐 website](https://rajewsky-lab.github.io/openst/latest) | [📜 preprint](https://www.biorxiv.org/content/10.1101/2023.12.22.572554v1) | [🐁 datasets](https://rajewsky-lab.github.io/openst/latest/examples/datasets/)

Open-ST is an open-source spatial transcriptomics method 
with efficient whole-transcriptome capture at sub-cellular resolution (0.6 μm) at low cost 
(<150 Euro library preparation per 12 mm²).

Open-ST requires standard lab equipment
and includes open-source software for seamless data processing and analysis.

- Discuss development on [GitHub](https://github.com/rajewsky-lab/openst).
- Read the [step-by-step protocol](https://rajewsky-lab.github.io/openst/latest/introduction/).
- Install via `pip install openst`.

## Tutorials and examples
We provide [various Open-ST datasets](https://rajewsky-lab.github.io/openst/examples/getting_started/) collected by the [Rajewsky lab @ MDC Berlin](https://www.mdc-berlin.de/n-rajewsky), from various tissues/organisms.
There, links to raw and processed data are available, as well as step-by-step guides for (pre)processing and downstream analysis. 

## Documentation
All the detail to the experimental and computational aspects of our method are available in [our website](https://rajewsky-lab.github.io/openst/).

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
   pip install mkdocs-material mkdocs-autorefs mknotebooks
   ```
3. Serve mkdocs with the following command:
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

## Code of Conduct
Everyone interacting in the Open-ST project's codebases, issue trackers, and discussion forums is expected to follow the [PSF Code of Conduct](https://www.python.org/psf/conduct/).

## License
The software tools of this project are under the GNU License - see the [LICENSE](LICENSE) file for details.
