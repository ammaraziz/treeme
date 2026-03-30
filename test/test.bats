#!/usr/bin/env bats

global_setup() {
    mkdir -p test/test-data/output/
    TREE_BASIC="test/test-data/test.10.nwk"
    TREE_BOOT="test/test-data/test.10.nwk"
    TREE_MED="test/test-data/test.100.nwk"
    TREE_NEXUS="test/test-data/test.10.nexus"
    META="test/test-data/test.10.tsv"
    OUTPUT_BASE="test/test-data/output/"
    OUTPUT_FAILED="test/test-data/output/fail.pdf"
}

global_teardown() {
    rm -f $OUTPUT_BASE/*.pdf
}

global_setup

@test "Test 1: missing argument - Failure" {
    run ./treeme.R -t $TREE_BASIC -o $OUTPUT_FAILED

    [[ $status -eq 1 ]]
    [[ "$output" =~ "Missing arguments: --meta" ]] # =~ is bats-core for "contains"
}

@test "Test 2: basic test - Success" {
    run ./treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test2_basic.pdf"

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test2_basic.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test2_basic.pdf" ]] 
}

@test "Test 3: basic tree with no labels - Success" {
    run ./treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test3_nolabs.pdf" --taxa-font-size 0

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test3_nolabs.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test3_nolabs.pdf" ]] 
}

@test "Test 4: color taxa - Success" {
    run ./treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test4_colortaxa.pdf" --color-by "epicluster"

    [[ "$status" -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test4_colortaxa.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test4_colortaxa.pdf" ]] 
}

@test "Test 5: color taxa - Failure" {
    run ./treeme.R -t $TREE_BASIC -m $META -o $OUTPUT_FAILED --color-by "NotReal"

    [[ "$status" -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]] # =~ is bats-core for "contains"
}

@test "Test 6: tree with tippoint - Success" {
    run ./treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test6_tippoint.pdf" --shape-by "vaccination_status"

    [[ "$status" -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test6_tippoint.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test6_tippoint.pdf" ]] 
}

@test "Test 7: tippoint missing shape-var - Failure" {
    run ./treeme.R -t $TREE_BASIC -m $META -o $OUTPUT_FAILED --shape-by "NotReal"

    [[ "$status" -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]]
}

@test "Test 8: colors and shape input - Success" {
    run ./treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test8_taxacolorshapes.pdf" \
    --color-by "epicluster" --shape-by "vaccination_status"

    [[ "$status" -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test8_taxacolorshapes.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test8_taxacolorshapes.pdf" ]]
}

@test "Test 9: basic tree with no labels, shapes + shape colors - Success" {
    run ./treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test9_shapecolor_nolabels.pdf" \
    --taxa-font-size 0 --color-by "epicluster" --shape-by "vaccination_status"

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test9_shapecolor_nolabels.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test9_shapecolor_nolabels.pdf" ]] 
}

@test "Test 11: no taxa, shapes - Success" {
    run ./treeme.R -t $TREE_BOOT -m $META -o "${OUTPUT_BASE}/test11_big_shapescolors_nolabels.pdf" \
    --taxa-font-size 0 --shape-by "collection_year"

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test11_big_shapescolors_nolabels.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test11_big_shapescolors_nolabels.pdf" ]] 
}

@test "Test 12: branch labels - Success" {
    run ./treeme.R -t $TREE_BOOT -m $META -o "${OUTPUT_BASE}/test12_bootstrap.pdf" --bootstrap

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test12_bootstrap.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test12_bootstrap.pdf" ]] 
}

@test "Test 13: med tree with labs, shapes, bootstrap" {
    run ./treeme.R -t $TREE_MED -m $META -o "${OUTPUT_BASE}/test13_medtree.pdf" \
    --color-by "epicluster" \
    --shape-by "vaccination_status"

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test12_bootstrap.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test12_bootstrap.pdf" ]] 
}

#@test "Test X: heatmap - Success" {
#    run ./treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test13_heatmap.pdf" \
#    --heatmap-value-range "7:8" --heatmap-fill-range "9:10"
#
#    [[ $status -eq 0 ]]
#    [[ -f "${OUTPUT_BASE}/test13_heatmap.pdf" ]]
#    [[ -s "${OUTPUT_BASE}/test13_heatmap.pdf" ]] 
#}

# global_teardown

