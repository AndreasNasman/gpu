const convertInputData = require("./input-data/convert-input-data");

const hrstart = process.hrtime();

convertInputData.execute();

const hrend = process.hrtime(hrstart);
console.info("Execution time (hr): %ds %dms\n", hrend[0], hrend[1] / 1000000);
