#!/usr/bin/env bats

setup() {
    TREE="test/test-data/basic.nwk"
    META="test/test-data/metadata.tsv"
    OUTPUT="test/test-data/output/basic.pdf"
    COLOR_FILE="test/test-data/colors.tsv"
    COLOR_VAR="month"
}

teardown() {
    rm -f $OUTPUT
}

@test "Test: basic test" {
    treeme.R -t $TREE -m $META -o $OUTPUT

    # assertions
    [[ $status -eq 0 ]]
    [[ -f "$OUTPUT" ]]
    [[ -s "$OUTPUT" ]] 
}

@test "Test: basic tree with no labels" {
    treeme.R -t $TREE -m $META -o $OUTPUT --font-size 0

    # assertions
    [[ $status -eq 0 ]]
    [[ -f "$OUTPUT" ]]
    [[ -s "$OUTPUT" ]] 
}

@test "Test: color taxa - failure" {
    treeme.R -t $TREE -m $META -o $OUTPUT --colors-file $COLOR_FILE --color-by-var $COLOR_VAR

    # assertions
    [[ $status -eq 0 ]]
    [[ -f "$OUTPUT" ]]
    [[ -s "$OUTPUT" ]]
}