exports.execute = () => {
  console.log("\x1b[33m%s\x1b[0m", "Conversion started.");

  const fs = require("fs");

  const arcminutesToRadians = arcminutes =>
    (Number(arcminutes) * Math.PI) / (60 * 180);

  const convertToRadiansJSON = (inputFile, outputName) => {
    const data = fs.readFileSync(
      `${__dirname}/original/${inputFile}.txt`,
      "utf-8"
    );

    const result = data
      .split("\n")
      .slice(1, -1)
      .map(row => {
        let [rightAscension, declination] = row.split(/\s+/);
        declination = arcminutesToRadians(declination);
        rightAscension = arcminutesToRadians(rightAscension);

        return {
          declination,
          rightAscension
        };
      });

    fs.writeFileSync(`${__dirname}/${outputName}.json`, JSON.stringify(result));
    console.log(`'${outputName}' was converted successfully!`);
  };

  convertToRadiansJSON("data_100k_arcmin", "real-galaxies");
  convertToRadiansJSON("flat_100k_arcmin", "random-galaxies");

  console.log("\x1b[33m%s\x1b[0m", "Conversion finished.\n");
};
