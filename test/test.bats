#!/usr/bin/env bats

setup() {
    TREE="test/test-data/basic.nwk"
    META="test/test-data/metadata.tsv"
    OUTPUT="test/test-data/output/basic.pdf"
    COLOR_FILE="test/test-data/colors.tsv"
    COLOR_VAR="month"
    SHAPES_FILE="test/test-data/shapes.tsv"
    SHAPE_VAR="state"
}

teardown() {
    rm -f $OUTPUT
}

@test "Test: missing argument - Failure" {
    run treeme.R -t $TREE -o $OUTPUT

    [[ $status -eq 1 ]]
    [[ "$output" =~ "Missing arguments: --meta" ]] # =~ is bats-core for "contains"
}

@test "Test: basic test - Success" {
    run treeme.R -t $TREE -m $META -o $OUTPUT

    # assertions
    [[ $status -eq 0 ]]
    [[ -f "$OUTPUT" ]]
    [[ -s "$OUTPUT" ]] 
}

@test "Test: basic tree with no labels - Success" {
    run treeme.R -t $TREE -m $META -o $OUTPUT --font-size 0

    # assertions
    [[ $status -eq 0 ]]
    [[ -f "$OUTPUT" ]]
    [[ -s "$OUTPUT" ]] 
}

@test "Test: color taxa - Success" {
    treeme.R -t $TREE -m $META -o $OUTPUT --colors-file $COLOR_FILE --color-by-var $COLOR_VAR

    [[ "$status" -eq 0 ]]
    [[ -f "$OUTPUT" ]]
    [[ -s "$OUTPUT" ]] 
}

@test "Test: color taxa - Failure" {
    run treeme.R -t $TREE -m $META -o $OUTPUT --colors-file $COLOR_FILE --color-by-var NotReal

    [[ "$status" -eq 1 ]]
    [[ "$output" =~ "CRITICAL" ]] # =~ is bats-core for "contains"
}

@test "Test: tippoint - Success" {
    treeme.R -t $TREE -m $META -o $OUTPUT --shapes-file $SHAPES_FILE --shape-by-var $SHAPE_VAR

    [[ "$status" -eq 0 ]]
    [[ -f "$OUTPUT" ]]
    [[ -s "$OUTPUT" ]] 
}