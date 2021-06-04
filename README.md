# mpn-aml-pairtree

Source code used in combination with Pairtree for studying the cancer evolution from myeloproliferative neoplasm (MPN) to acute myeloid leukemia (AML).


## Requirements
Python 3+


## Installation

Clone repository
```
git clone https://github.com/ethanumn/mpn-aml-pairtree.git && cd mpn-aml-pairtree
```

Setup virtual environment

```
virtualenv env
echo "" >> env/bin/activate
echo "# Environment Variables" >> env/bin/activate
echo "export DATA_DIR=$PWD/data" >> env/bin/activate
echo "export UTILS_DIR=$PWD/utils" >> env/bin/activate
source env/bin/activate
```

Install dependencies

```
pip install -r requirements.txt
```
