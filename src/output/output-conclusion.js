exports.execute = () => {
  console.log("\x1b[33m%s\x1b[0m", "Conlusion started.");

  const fs = require("fs");

  const histogramDD = require("../histograms/DD.json");
  const histogramDR = require("../histograms/DR.json");
  const histogramRR = require("../histograms/RR.json");

  const result = [];
  for (let i = 0; i < histogramDD.length; i += 1) {
    result.push(
      (histogramDD[i] - 2 * histogramDR[i] + histogramRR[i]) / histogramRR[i]
    );
  }

  fs.writeFileSync(`${__dirname}/result.json`, JSON.stringify(result));

  console.log("\x1b[33m%s\x1b[0m", "Conlusion finished.");
};
