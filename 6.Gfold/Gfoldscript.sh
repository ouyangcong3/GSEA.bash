 ./gfold/gfold diff -s1 HCT116_NC1.cnt -s2 HCT116_si1.cnt -o rep1.diff
./gfold/gfold diff -s1 HCT116_NC2.cnt -s2 HCT116_si2.cnt -o rep2.diff
cat rep1.diff | tail -n +11 | awk '{print $1"\t"$5}' >siPHLDA1_rep1_diff.rnk
cat rep2.diff | tail -n +11 | awk '{print $1"\t"$5}' >siPHLDA1_rep2_diff.rnk