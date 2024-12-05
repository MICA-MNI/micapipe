snakemake --dag > dag.dot 
dot -Tsvg dag.dot -o dag.svg    
rm dag.dot