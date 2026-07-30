#!/bin/bash
exp_heatmap \
    benchmark \
    -s 40000000 \
    -e 40500000 \
    -o chr22_bm_1 \
    --report chr22_bm_1_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr22_p3_37.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/testing_changes/genotypes.panel

exp_heatmap \
    benchmark \
    -s 16060513 \
    -e 51210268 \
    -o chr22_bm_2 \
    --report chr22_bm_2_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr22_p3_37.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/testing_changes/genotypes.panel

exp_heatmap \
    benchmark \
    -s 40000000 \
    -e 40500000 \
    -r 5 \
    -o chr22_bm_3 \
    --report chr22_bm_3_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr22_p3_37.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/testing_changes/genotypes.panel

exp_heatmap \
    benchmark \
    -s 16060513 \
    -e 51210268 \
    -r 5 \
    -o chr22_bm_4 \
    --report chr22_bm_4_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr22_p3_37.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/testing_changes/genotypes.panel

exp_heatmap \
    benchmark \
    -s 40000000 \
    -e 40500000 \
    -r 5 \
    -w 1 \
    -o chr22_bm_5 \
    --report chr22_bm_5_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr22_p3_37.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/testing_changes/genotypes.panel

exp_heatmap \
    benchmark \
    -s 16060513 \
    -e 51210268 \
    -r 5 \
    -w 1 \
    -o chr22_bm_6 \
    --report chr22_bm_6_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr22_p3_37.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/testing_changes/genotypes.panel
