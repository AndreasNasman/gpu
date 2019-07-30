const hrstart = process.hrtime();

const convertInputData = require("./input-data/convert-input-data");
convertInputData.execute();

const buildHistograms = require("./histograms/build-histograms");
buildHistograms.execute();

const outputConclusion = require("./output/output-conclusion");
outputConclusion.execute();

const hrend = process.hrtime(hrstart);
console.info("\nExecution time (hr): %ds %dms", hrend[0], hrend[1] / 1000000);
