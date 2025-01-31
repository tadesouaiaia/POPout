# POPout
The POPOut Test 


# To run POPout in Python: 


./POPout.py --tailSize 1 -o test data/toy_bp.pop data/toy_vitD.pop data/toy_hg.pop data/toy_bmi.pop data/toy_height.pop

# To create a POP file using Python: 


./POPout.py --makePOPfile data/TEST.prs data/TEST.pheno -o my_file




# To run POPout in R (on a single file): 

# Rscript --vanilla POPout.R data/toy_bp.pop 

# To run POPout in R on separate files: 

# Rscript --vanilla POPout.R data/TEST.prs data/TEST.pheno

