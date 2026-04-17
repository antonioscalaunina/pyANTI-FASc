docker run --rm -it `
  -e MPLCONFIGDIR=/tmp/matplotlib `
  -v "${PWD}:/app" `
  pyantifasc
