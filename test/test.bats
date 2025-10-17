#!/usr/bin/env bats

global_setup() {
    TREE_BASIC="test/test-data/basic.nwk"
    TREE_BOOT="test/test-data/boot.nwk"
    TREE_NEXUS="test/test-data/test.nexus"
    META="test/test-data/metadata.tsv"
    META_BOOT="test/test-data/metadata.boot.tsv"
    OUTPUT_BASE="test/test-data/output/"
    OUTPUT_FAILED="test/test-data/output/fail.pdf"
}

global_teardown() {
    rm -f $OUTPUT_BASE/*.pdf
}

global_setup

@test "Test: missing argument - Failure" {
    run treeme.R -t $TREE_BASIC -o $OUTPUT_FAILED

    [[ $status -eq 1 ]]
    [[ "$output" =~ "Missing arguments: --meta" ]] # =~ is bats-core for "contains"
}

@test "Test: basic test - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/basic.pdf"

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/basic.pdf" ]]
    [[ -s "${OUTPUT_BASE}/basic.pdf" ]] 
}

@test "Test: basic tree with no labels - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/notext.pdf" --font-size 0

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/notext.pdf" ]]
    [[ -s "${OUTPUT_BASE}/notext.pdf" ]] 
}

@test "Test: color taxa - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/colortaxa.pdf" -c "month"

    [[ "$status" -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/colortaxa.pdf" ]]
    [[ -s "${OUTPUT_BASE}/colortaxa.pdf" ]] 
}

@test "Test: color taxa - Failure" {
    run treeme.R -t $TREE_BASIC -m $META -o $OUTPUT_FAILED -c "NotReal"

    [[ "$status" -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]] # =~ is bats-core for "contains"
}

@test "Test: tippoint - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o $OUTPUT_BASE -s "state"

    [[ "$status" -eq 0 ]]
    [[ -f "$OUTPUT_BASE" ]]
    [[ -s "$OUTPUT_BASE" ]] 
}

@test "Test: tippoint missing shape-var - Failure" {
    run treeme.R -t $TREE_BASIC -m $META -o $OUTPUT_FAILED -s "NotReal"

    [[ "$status" -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]]
}

@test "Test: colors and shape input - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/colorshape.pdf" \
    -c "month" -s "state"

    [[ "$status" -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/colorshape.pdf" ]]
    [[ -s "${OUTPUT_BASE}/colorshape.pdf" ]]
}

@test "Test: basic tree with no labels, has colors and shapes - Success" {
    run treeme.R -t $TREE_BASIC -m $META -o "${OUTPUT_BASE}/shapecolor_notext.pdf" \
    --font-size 0 -c "month" -s "state"

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/shapecolor_notext.pdf" ]]
    [[ -s "${OUTPUT_BASE}/shapecolor_notext.pdf" ]] 
}

@test "Test: big tree; no taxa with shapes - Failed" {
    run treeme.R -t $TREE_BOOT -m $META_BOOT -o "${OUTPUT_BASE}/boot_shapecolor.pdf" \
    --font-size 0 -c "month" -s "state"

    [[ $status -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]]

}

@test "Test: big tree; no taxa with shapes - Success" {
    run treeme.R -t $TREE_BOOT -m $META_BOOT -o "${OUTPUT_BASE}/boot_shapecolor.pdf" \
    --font-size 0 -s "state"

    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/boot_shapecolor.pdf" ]]
    [[ -s "${OUTPUT_BASE}/boot_shapecolor.pdf" ]] 
}

#global_teardown