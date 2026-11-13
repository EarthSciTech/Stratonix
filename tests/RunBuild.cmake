# Run `stratonix build` into a fresh folder and fail unless it exits 0 (all checks pass).
file(REMOVE_RECURSE ${OUT})
execute_process(COMMAND ${EXE} build ${DESC} -o ${OUT} --no-figures RESULT_VARIABLE rc OUTPUT_VARIABLE so ERROR_VARIABLE se)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "stratonix build failed (${rc}):\n${so}\n${se}")
endif()
