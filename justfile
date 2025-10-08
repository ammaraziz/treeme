# Test treeme.r CLI tool
treeme-test:
    treeme.r -t test-data/test.nwk -o output.pdf -m test-data/metadata.tsv

    # if test -f output.pdf; then
    #     echo "PASS: PDF created successfully"
    #     exit 0
    # else
    #     echo "FAIL: PDF not created"
    #     exit 1
    # fi
    echo building