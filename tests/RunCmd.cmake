# Run a stratonix command into a fresh output folder; fail unless it exits 0.
# -DEXE=... -DOUT=<folder> -DARGS="cmd;arg1;arg2..." (the output folder is appended as -o OUT)
file(REMOVE_RECURSE ${OUT})
file(REMOVE ${OUT}.zip)
execute_process(COMMAND ${EXE} ${ARGS} -o ${OUT} RESULT_VARIABLE rc OUTPUT_VARIABLE so ERROR_VARIABLE se)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "stratonix ${ARGS} failed (${rc}):\n${so}\n${se}")
endif()
