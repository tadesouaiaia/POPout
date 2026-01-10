# POPout
The POPOut Test 


### PYTHON3 Commands ###

# 1) To run POPout on a Single POP File

./POPout.py --tailSize 1 -o test data/toy_bp.pop --name "Blood Pressure"  


# 2) To run POPout on a Multiple POP Files

./POPout.py --tailSize 1 -o test data/toy_bp.pop data/toy_vitD.pop data/toy_hg.pop data/toy_bmi.pop data/toy_height.pop


# 3a) To Run POPout on a Multiple POP Files and Select Names and Colors (for forest plot) 

./POPout.py -o test data/toy_bp.pop ../../data/toy_height.pop --names "Blood Pressure" "Height" --colors red orange


# 3b) Alternatively, set names/colors using a config file: 

./POPout.py -o test data/*.pop --config ../data/toy.config


# 4) To create a POP file using Python: 

./POPout.py --makePOPfile data/TEST.prs data/TEST.pheno -o my_file


#######################################################




### R Commands ###

# To run POPout in R (on a single file): 

# Rscript --vanilla POPout.R data/toy_bp.pop 

# To run POPout in R on separate files: 

# Rscript --vanilla POPout.R data/TEST.prs data/TEST.pheno

