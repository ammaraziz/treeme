# Test treeme.r CLI tool
treeme-test:
    treeme.r -t test-data/test.nwk -o test-data/output.pdf -m test-data/metadata.tsv --paper-size A5l --font-size 0

    # if test -f output.pdf; then
    #     echo "PASS: PDF created successfully"
    #     exit 0
    # else
    #     echo "FAIL: PDF not created"
    #     exit 1
    # fi
    echo building