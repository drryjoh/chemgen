# Chemgen

A brief description of your project goes here.

## Table of Contents

- [Installation](#installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [License](#license)

## Installation

To set up the project on your local machine, please follow these steps:

1. Clone this repository:

   ```bash
   git clone https://github.com/drryjoh/chemgen.git
   cd chemgen
   ```

2. Install the required dependencies using Python 3:

We recommend installing a unique python environment using `python3.11` as that is the latest python version that supports cantera. To create a unique environmnet:

```terminal
mkdir ~/python_environments/
cd ~/python_environments/
python3.11 -m venv chemgen
source ~/python_environments/chemgen/bin/activate
```

you can add to you bashrc/zshrc file
```
alias chemgen="source ~/python_environments/chemgen/bin/activate"
```

which will allow you to access the source from command line

```terminal
chemgen
```

Once you have your unique python environment, you can install the requirements via

   ```bash
   python3 -m pip install -r requirements.txt
   ```

   If you don’t have a `requirements.txt` file, you can manually install the dependencies as follows:

   ```bash
   python3 -m pip install cantera pyyaml
   ```

## Usage

The best way to learn to use ChemGen is through the tutorials located [here](tutorial/README.md).

## License


