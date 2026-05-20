#!/bin/bash
LOG_FILE="/scratch/groups/wangbo/Multiome_test/work/df/3af2520296088305d06cba0b412151/e18_test_multiome.chromap.log"
echo "--- CHROMAP LOG ---" > chromap_error.txt
if [ -f "$LOG_FILE" ]; then cat "$LOG_FILE" >> chromap_error.txt; else echo "NOT FOUND" >> chromap_error.txt; fi
echo "--- COMMAND ERR ---" >> chromap_error.txt
cat /scratch/groups/wangbo/Multiome_test/work/df/3af2520296088305d06cba0b412151/.command.err >> chromap_error.txt 2>/dev/null || true
echo "--- COMMAND OUT ---" >> chromap_error.txt
cat /scratch/groups/wangbo/Multiome_test/work/df/3af2520296088305d06cba0b412151/.command.out >> chromap_error.txt 2>/dev/null || true
