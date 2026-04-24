$MODE = $args[0]

if ($MODE -eq "notebook") {
    docker run --rm -it `
        -e HOME=/tmp `
        -e JUPYTER_RUNTIME_DIR=/tmp/jupyter-runtime `
        -e JUPYTER_DATA_DIR=/tmp/jupyter-data `
        -e JUPYTER_CONFIG_DIR=/tmp/jupyter-config `
        -e MPLCONFIGDIR=/tmp/matplotlib `
        -p 8888:8888 `
        -v "${PWD}:/app" `
        pyantifasc `
        conda run --no-capture-output -n antifasc jupyter lab `
        --ip=0.0.0.0 --port=8888 --no-browser
}
else {
    docker run --rm -it `
        -e MPLCONFIGDIR=/tmp/matplotlib `
        -v "${PWD}:/app" `
        pyantifasc `
        conda run --no-capture-output -n antifasc python $args
}
