# This will be for analysis of data gathered from XRD

using DataFrames, CSV

"""
This function reads the CSV file located at `fileloc`.
The reader should know whether the file location given is absolute or relative.
"""
function readxrd(fileloc::String)
  return CSV.read(fileloc, header=["theta", "counts"])
end