#!/bin/bash
exp_heatmap \
    benchmark \
    -s 136108646 \
    -e 137108646 \
    -o chr2_bm_1 \
    -c \
    --report chr2_bm_1_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr2_38.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/genes_testing/LCT_mytake_38/30xfinal.panel

exp_heatmap \
    benchmark \
    -s 47920990 \
    -e 48920990 \
    -o chr2_bm_2 \
    -c \
    --report chr2_bm_2_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr2_38.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/genes_testing/LCT_mytake_38/30xfinal.panel

exp_heatmap \
    benchmark \
    -s 136108646 \
    -e 137108646 \
    -r 5 \
    -w 1 \
    -o chr2_bm_3 \
    -c \
    --report chr2_bm_3_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr2_38.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/genes_testing/LCT_mytake_38/30xfinal.panel

exp_heatmap \
    benchmark \
    -s 47920990 \
    -e 48920990 \
    -r 5 \
    -w 1 \
    -o chr2_bm_4 \
    -c \
    --report chr2_bm_4_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr2_38.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/genes_testing/LCT_mytake_38/30xfinal.panel

exp_heatmap \
    benchmark \
    -s 40000000 \
    -e 137000000 \
    -o chr2_bm_5 \
    -c \
    --report chr2_bm_5_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr2_38.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/genes_testing/LCT_mytake_38/30xfinal.panel

exp_heatmap \
    benchmark \
    -s 40000000 \
    -e 137000000 \
    -r 5 \
    -w 1 \
    -o chr2_bm_6 \
    -c \
    --report chr2_bm_6_report \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/vcf_files/chr2_38.vcf \
    /home/adam/Desktop/PhD/Cursor_projects/exp_heatmap/tool_testing/genes_testing/LCT_mytake_38/30xfinal.panel