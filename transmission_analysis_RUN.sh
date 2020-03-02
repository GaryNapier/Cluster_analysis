#!/bin/bash


# ----------TESTING ------------------
# python_scripts/transmission_analysis.py -l 4.1 \
# -dm lin_4_results.filt.dist \
#  -mds lin_4_results.filt.pca.mds \
#   -lf lin_4.1_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt

# python_scripts/transmission_analysis.py -l 1.2 \
# -dm lin_1_results.filt.dist \
#  -mds lin_1_results.filt.pca.mds \
#   -lf lin_1.2_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# ----------TESTING ------------------

COUNTER=1
# for lineage in $"1" $"1.1" $"1.2" $"2" $"3" $"4" $"4.1" $"4.2" $"4.3" $"4.4" $"4.5" $"4.6" $"4.7" $"4.8" $"4.9" $"4.10" $"4.11" $"5" $"6" $"7" $"BOV"
for lineage in $"1" $"2" $"3" $"4" $"5" $"6" $"7" $"BOV"
do
# for clust in "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L"
for sublin in $"1" $"1.1" $"1.1.1" $"1.1.1.1" $"1.1.2" $"1.1.3" \
$"1.2" $"1.2.1" $"1.2.2" \
 $"2" $"2.1" $"2.2.1" $"2.2.2" \
 $"3" $"3.1" $"3.1.1" $"3.1.2" $"3.1.2.1" $"3.1.2.2" \
 $"4" $"4.1" $"4.1.1" $"4.1.1.1" $"4.1.1.2" $"4.1.1.3" $"4.1.2" $"4.1.2.1" \
 $"4.2" $"4.2.1" $"4.2.2" $"4.2.2.1" \
 $"4.3" $"4.3.1" $"4.3.2" $"4.3.2.1" $"4.3.3" $"4.3.4" $"4.3.4.1" $"4.3.4.2" $"4.3.4.2.1" \
 $"4.4" $"4.4.1" $"4.4.1.1" $"4.4.1.2" $"4.4.2" \
 $"4.5" \
 $"4.6" $"4.6.1" $"4.6.1.1" $"4.6.1.2" $"4.6.2" $"4.6.2.1" $"4.6.2.2" \
 $"4.7" \
 $"4.8" \
 $"4.9" \
  $"4.10" $"4.11" \
	$"5" $"6" $"7" $"BOV"
do
# if [ -f lin_${lineage}/txt/lin_${sublin}_true_samps.txt ]; then
echo ${lineage}
echo ${sublin}
echo " Loop number: $COUNTER "
COUNTER=$[$COUNTER +1]
python_scripts/transmission_analysis.py -l ${sublin} \
-dm lin_${lineage}_results.filt.dist \
 -mds lin_${lineage}_results.filt.pca.mds \
  -lf lin_${sublin}_true_samps.txt \
   -ff revised_sublins_fst_table.txt \
   -trans_table_out all_lins/txt/transmission_table.txt  \
   -trans_samps_out all_lins/txt/transmission_samples_list.txt
# fi
done
done

# for lineage in $"1"
# do
# for sublin in $"1"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"1.1"
# do
# for sublin in $"1.1" $"1.1.1" $"1.1.1.1" $"1.1.2" $"1.1.3"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"1.2"
# do
# for sublin in $"1.2" $"1.2.1" $"1.2.2"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"2"
# do
# for sublin in  $"2" $"2.1" $"2.2.1" $"2.2.2"
# do
#   python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"3"
# do
# for sublin in  $"3" $"3.1" $"3.1.1" $"3.1.2" $"3.1.2.1" $"3.1.2.2"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4"
# do
# for sublin in  $"4"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4.1"
# do
# for sublin in  $"4.1" $"4.1.1" $"4.1.1.1" $"4.1.1.2" $"4.1.1.3" $"4.1.2" $"4.1.2.1"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4.2"
# do
# for sublin in  $"4.2" $"4.2.1" $"4.2.2" $"4.2.2.1"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4.3"
# do
# for sublin in $"4.3" $"4.3.1" $"4.3.2" $"4.3.2.1" $"4.3.3" $"4.3.4" $"4.3.4.1" $"4.3.4.2" $"4.3.4.2.1"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4.4"
# do
# for sublin in $"4.4" $"4.4.1" $"4.4.1.1" $"4.4.1.2" $"4.4.2"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4.5"
# do
# for sublin in $"4.5"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4.6"
# do
# for sublin in $"4.6" $"4.6.1" $"4.6.1.1" $"4.6.1.2" $"4.6.2" $"4.6.2.1" $"4.6.2.2"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4.7"
# do
# for sublin in $"4.7"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4.8"
# do
# for sublin in $"4.8"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4.9"
# do
# for sublin in $"4.9"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4"
# do
# for sublin in $"4.10"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"4"
# do
# for sublin in $"4.11"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"5"
# do
# for sublin in $"5"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"6"
# do
# for sublin in $"6"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"7"
# do
# for sublin in $"7"
# do
# python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
#
# for lineage in $"BOV"
# do
# for sublin in $"BOV"
# do
#   python_scripts/transmission_analysis.py -l ${sublin} \
# -dm lin_${lineage}_results.filt.dist \
#  -mds lin_${lineage}_results.filt.pca.mds \
#   -tf lin_${lineage}.iqtree.treefile \
#   -lf lin_${sublin}_true_samps.txt \
#    -ff revised_sublins_fst_table.txt \
#    -trans_table_out all_lins/txt/transmission_table.txt  \
#    -trans_samps_out all_lins/txt/transmission_samples_list.txt
# done
# done
