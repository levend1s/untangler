# untangler
Assembler for DNA, aligner for RNA against DNA, untangler for RNA subclusters

# SETUP
```
# set up and enter virtual environment
$ python3 -m venv .venv
$ source .venv/bin/activate
$ pip install -r requirements.txt

or

$ python3 -m virtualenv env
$ source env/bin/activate
$ python -m pip install -r requirements.txt

# run untangler

# exit the virtual environment
$ deactivate
```

# Usage

```
$ ID=PF3D7_1007900; python3 untangler.py cluster_transcripts --ids ${ID} -i config.txt -o ./cluster_transcripts_${ID}_28 --cluster_cols introns,m6A --predict introns -t 0.6 -f --generate_igv_reports 
```