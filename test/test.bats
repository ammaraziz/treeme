#!/usr/bin/env bats

global_setup() {
    TREE="test/test-data/basic.nwk"
    META="test/test-data/metadata.tsv"
    OUTPUT_BASE="test/test-data/output/"
    OUTPUT_FAILED="test/test-data/output/fail.pdf"
    COLOR_FILE="test/test-data/colors.tsv"
    COLOR_VAR="month"
    SHAPES_FILE="test/test-data/shapes.tsv"
    SHAPE_VAR="state"
}

global_teardown() {
    rm -f $OUTPUT_BASE/*.pdf
}

global_setup

@test "Test: missing argument - Failure" {
    run treeme.R -t $TREE -o $OUTPUT_FAILED

    [[ $status -eq 1 ]]
    [[ "$output" =~ "Missing arguments: --meta" ]] # =~ is bats-core for "contains"
}

@test "Test: basic test - Success" {
    run treeme.R -t $TREE -m $META -o "${OUTPUT_BASE}/basic.pdf"

    # assertions
    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/basic.pdf" ]]
    [[ -s "${OUTPUT_BASE}/basic.pdf" ]] 
}

@test "Test: basic tree with no labels - Success" {
    run treeme.R -t $TREE -m $META -o "${OUTPUT_BASE}/notext.pdf" --font-size 0

    # assertions
    [[ $status -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/notext.pdf" ]]
    [[ -s "${OUTPUT_BASE}/notext.pdf" ]] 
}

@test "Test: color taxa - Success" {
    run treeme.R -t $TREE -m $META -o "${OUTPUT_BASE}/colortaxa.pdf" \
    --colors-file $COLOR_FILE --color-by-var $COLOR_VAR

    [[ "$status" -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/colortaxa.pdf" ]]
    [[ -s "${OUTPUT_BASE}/colortaxa.pdf" ]] 
}

@test "Test: color taxa - Failure" {
    run treeme.R -t $TREE -m $META -o $OUTPUT_FAILED \
    --colors-file $COLOR_FILE --color-by-var "NotReal"

    [[ "$status" -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]] # =~ is bats-core for "contains"
}

@test "Test: tippoint - Success" {
    run treeme.R -t $TREE -m $META -o $OUTPUT_BASE \
    --shapes-file $SHAPES_FILE --shape-by-var $SHAPE_VAR

    [[ "$status" -eq 0 ]]
    [[ -f "$OUTPUT_BASE" ]]
    [[ -s "$OUTPUT_BASE" ]] 
}

@test "Test: tippoint missing shape-file - Failure" {
    run treeme.R -t $TREE -m $META -o $OUTPUT_FAILED \
    --shapes-file $SHAPES_FILE

    [[ "$status" -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]]
}

@test "Test: tippoint missing shape-var - Failure" {
    run treeme.R -t $TREE -m $META -o $OUTPUT_FAILED \
    --shape-by-var $SHAPE_VAR

    [[ "$status" -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]]
}

@test "Test: colors and shape input - Success" {
    run treeme.R -t $TREE -m $META -o "${OUTPUT_BASE}/colorshape.pdf" \
    --shapes-file $SHAPES_FILE --shape-by-var $SHAPE_VAR \
    --colors-file $COLOR_FILE --color-by-var $COLOR_VAR

    [[ "$status" -eq 0 ]]
    [[ -f "${OUTPUT_BASE}/colorshape.pdf" ]]
    [[ -s "${OUTPUT_BASE}/colorshape.pdf" ]] 
}

global_teardown