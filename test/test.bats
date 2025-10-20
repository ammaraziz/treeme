#!/usr/bin/env bats

global_setup() {
    TREE_BASIC="test/test-data/test.nwk"
    TREE_BOOT="test/test-data/boot.nwk"
    TREE_NEXUS="test/test-data/test.nexus"
    TREE_NHX="test/test-data/test.nhx"
    META="test/test-data/metadata.tsv"
    META_BOOT="test/test-data/metadata.boot.tsv"
    OUTPUT_BASE="test/test-data/output/"
    OUTPUT_FAILED="test/test-data/output/fail.pdf"
}

global_teardown() {
    rm -f $OUTPUT_BASE/*.pdf
}

global_setup

@test "Test 1: missing argument - Failure" {
    run treeme.R -t $TREE_BASIC -o $OUTPUT_FAILED

    [[ $status -eq 1 ]]
    [[ "$output" =~ "Missing arguments: --meta" ]] # =~ is bats-core for "contains"
}

@test "Test 2: basic test - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test2_basic.pdf"

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test2_basic.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test2_basic.pdf" ]] 
}

@test "Test 3: basic tree with no labels - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test3_nolabs.pdf" --taxa-font-size 0

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test3_nolabs.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test3_nolabs.pdf" ]] 
}

@test "Test 4: color taxa - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test4_colortaxa.pdf" -c "month"

    [[ "$status" -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test4_colortaxa.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test4_colortaxa.pdf" ]] 
}

@test "Test 5: color taxa - Failure" {
    run treeme.R -t $TREE_BASIC -m $META -o $OUTPUT_FAILED -c "NotReal"

    [[ "$status" -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]] # =~ is bats-core for "contains"
}

@test "Test 6: tree with tippoint - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test6_tippoint.pdf" \
    -s "state"

    [[ "$status" -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test6_tippoint.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test6_tippoint.pdf" ]] 
}

@test "Test 7: tippoint missing shape-var - Failure" {
    run treeme.R -t $TREE_BASIC -m $META -o $OUTPUT_FAILED -s "NotReal"

    [[ "$status" -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]]
}

@test "Test 8: colors and shape input - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test8_taxacolorshapes.pdf" \
    -c "month" -s "state"

    [[ "$status" -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test8_taxacolorshapes.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test8_taxacolorshapes.pdf" ]]
}

@test "Test 9: basic tree with no labels, has colors and shapes - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/test9_shapecolor_nolabels.pdf" \
    --taxa-font-size 0 -c "month" -s "state"

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test9_shapecolor_nolabels.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test9_shapecolor_nolabels.pdf" ]] 
}

@test "Test 10: big tree; no taxa with shapes - Failed" {
    run treeme.R -t $TREE_BOOT -m $META_BOOT -o "${OUTPUT_BASE}/failed.pdf" \
    --taxa-font-size 0 -c "month" -s "state"

    [[ $status -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]]
}

@test "Test 11: big tree; no taxa with shapes - Success" {
    run treeme.R -t $TREE_BOOT -m $META_BOOT -o "${OUTPUT_BASE}/test11_big_shapescolors_nolabels.pdf" \
    --taxa-font-size 0 -s "state"

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test11_big_shapescolors_nolabels.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test11_big_shapescolors_nolabels.pdf" ]] 
}

@test "Test 12: big tree; branch labels - Success" {
    run treeme.R -t $TREE_BOOT -m $META_BOOT -o "${OUTPUT_BASE}/test12_bootstrap.pdf" \
    --bootstrap

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/test12_bootstrap.pdf" ]]
    [[ -s "${OUTPUT_BASE}/test12_bootstrap.pdf" ]] 
}

#global_teardown