// Builds Chapter 1 manuscript (Bentham Books format) as a .docx file.
// Usage: node build_chapter1.js
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, ImageRun, Table, TableRow, TableCell,
  AlignmentType, HeadingLevel, WidthType, ShadingType, BorderStyle, Footer,
  PageNumber, LevelFormat,
} = require("docx");

const FONT = "Times New Roman";
const FIG = path.join(__dirname, "figures");

// ---------------------------------------------------------------------------
// References (Bentham style). Numbered automatically in order of first citation.
// ---------------------------------------------------------------------------
const REFS = {
  godfray: "Godfray HCJ, Beddington JR, Crute IR, et al. Food security: The challenge of feeding 9 billion people. Science 2010; 327(5967): 812-8. [http://dx.doi.org/10.1126/science.1185383]",
  un2022: "United Nations, Department of Economic and Social Affairs, Population Division. World Population Prospects 2022: Summary of Results. New York: United Nations 2022.",
  fao2009: "FAO. How to Feed the World in 2050. Rome: Food and Agriculture Organization of the United Nations 2009.",
  tilman: "Tilman D, Balzer C, Hill J, Befort BL. Global food demand and the sustainable intensification of agriculture. Proc Natl Acad Sci USA 2011; 108(50): 20260-4. [http://dx.doi.org/10.1073/pnas.1116437108]",
  liakos: "Liakos KG, Busato P, Moshou D, Pearson S, Bochtis D. Machine learning in agriculture: A review. Sensors (Basel) 2018; 18(8): 2674. [http://dx.doi.org/10.3390/s18082674]",
  kamilaris: "Kamilaris A, Prenafeta-Boldú FX. Deep learning in agriculture: A survey. Comput Electron Agric 2018; 147: 70-90. [http://dx.doi.org/10.1016/j.compag.2018.02.016]",
  jha: "Jha K, Doshi A, Patel P, Shah M. A comprehensive review on automation in agriculture using artificial intelligence. Artif Intell Agric 2019; 2: 1-12. [http://dx.doi.org/10.1016/j.aiia.2019.05.004]",
  lowder: "Lowder SK, Skoet J, Raney T. The number, size, and distribution of farms, smallholder farms, and family farms worldwide. World Dev 2016; 87: 16-29. [http://dx.doi.org/10.1016/j.worlddev.2015.10.041]",
  gebbers: "Gebbers R, Adamchuk VI. Precision agriculture and food security. Science 2010; 327(5967): 828-31. [http://dx.doi.org/10.1126/science.1183899]",
  wolfert: "Wolfert S, Ge L, Verdouw C, Bogaardt MJ. Big data in smart farming – A review. Agric Syst 2017; 153: 69-80. [http://dx.doi.org/10.1016/j.agsy.2017.01.023]",
  zhangN: "Zhang N, Wang M, Wang N. Precision agriculture—a worldwide overview. Comput Electron Agric 2002; 36(2-3): 113-32. [http://dx.doi.org/10.1016/S0168-1699(02)00096-0]",
  mulla: "Mulla DJ. Twenty five years of remote sensing in precision agriculture: Key advances and remaining knowledge gaps. Biosyst Eng 2013; 114(4): 358-71. [http://dx.doi.org/10.1016/j.biosystemseng.2012.08.009]",
  klerkx: "Klerkx L, Jakku E, Labarthe P. A review of social science on digital agriculture, smart farming and agriculture 4.0: New contributions and a future research agenda. NJAS Wagening J Life Sci 2019; 90-91: 100315. [http://dx.doi.org/10.1016/j.njas.2019.100315]",
  basso: "Basso B, Antle J. Digital agriculture to design sustainable agricultural systems. Nat Sustain 2020; 3(4): 254-6. [http://dx.doi.org/10.1038/s41893-020-0510-0]",
  rose: "Rose DC, Chilvers J. Agriculture 4.0: Broadening responsible innovation in an era of smart farming. Front Sustain Food Syst 2018; 2: 87. [http://dx.doi.org/10.3389/fsufs.2018.00087]",
  verdouw: "Verdouw C, Tekinerdogan B, Beulens A, Wolfert S. Digital twins in smart farming. Agric Syst 2021; 189: 103046. [http://dx.doi.org/10.1016/j.agsy.2020.103046]",
  russell: "Russell S, Norvig P. Artificial Intelligence: A Modern Approach. 4th ed. Hoboken, NJ: Pearson 2020.",
  mccarthy: "McCarthy J, Minsky ML, Rochester N, Shannon CE. A proposal for the Dartmouth summer research project on artificial intelligence, August 31, 1955. AI Mag 2006; 27(4): 12-4.",
  mckinion: "McKinion JM, Lemmon HE. Expert systems for agriculture. Comput Electron Agric 1985; 1(1): 31-40. [http://dx.doi.org/10.1016/0168-1699(85)90004-3]",
  zadeh: "Zadeh LA. Fuzzy sets. Inf Control 1965; 8(3): 338-53. [http://dx.doi.org/10.1016/S0019-9958(65)90241-X]",
  cortes: "Cortes C, Vapnik V. Support-vector networks. Mach Learn 1995; 20(3): 273-97. [http://dx.doi.org/10.1007/BF00994018]",
  breiman: "Breiman L. Random forests. Mach Learn 2001; 45(1): 5-32. [http://dx.doi.org/10.1023/A:1010933404324]",
  sutton: "Sutton RS, Barto AG. Reinforcement Learning: An Introduction. 2nd ed. Cambridge, MA: MIT Press 2018.",
  lecun: "LeCun Y, Bengio Y, Hinton G. Deep learning. Nature 2015; 521(7553): 436-44. [http://dx.doi.org/10.1038/nature14539]",
  krizhevsky: "Krizhevsky A, Sutskever I, Hinton GE. ImageNet classification with deep convolutional neural networks. Adv Neural Inf Process Syst 2012; 25: 1097-105.",
  he2016: "He K, Zhang X, Ren S, Sun J. Deep residual learning for image recognition. Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR). Las Vegas, NV, USA. 2016; pp. 770-8. [http://dx.doi.org/10.1109/CVPR.2016.90]",
  hochreiter: "Hochreiter S, Schmidhuber J. Long short-term memory. Neural Comput 1997; 9(8): 1735-80. [http://dx.doi.org/10.1162/neco.1997.9.8.1735]",
  vaswani: "Vaswani A, Shazeer N, Parmar N, et al. Attention is all you need. Adv Neural Inf Process Syst 2017; 30: 5998-6008.",
  redmon: "Redmon J, Divvala S, Girshick R, Farhadi A. You only look once: Unified, real-time object detection. Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR). Las Vegas, NV, USA. 2016; pp. 779-88. [http://dx.doi.org/10.1109/CVPR.2016.91]",
  li2024: "Li J, Xu M, Xiang L, Chen D, Zhuang W, Yin X, et al. Foundation models in smart agriculture: Basics, opportunities, and challenges. Comput Electron Agric 2024; 222: 109032. [http://dx.doi.org/10.1016/j.compag.2024.109032]",
  kirillov: "Kirillov A, Mintun E, Ravi N, et al. Segment anything. Proceedings of the IEEE/CVF International Conference on Computer Vision (ICCV). Paris, France. 2023; pp. 4015-26.",
  tzounis: "Tzounis A, Katsoulas N, Bartzanas T, Kittas C. Internet of Things in agriculture, recent advances and future challenges. Biosyst Eng 2017; 164: 31-48. [http://dx.doi.org/10.1016/j.biosystemseng.2017.09.007]",
  elijah: "Elijah O, Rahman TA, Orikumhi I, Leow CY, Hindia MN. An overview of Internet of Things (IoT) and data analytics in agriculture: Benefits and challenges. IEEE Internet Things J 2018; 5(5): 3758-73. [http://dx.doi.org/10.1109/JIOT.2018.2844296]",
  zhangK: "Zhang C, Kovacs JM. The application of small unmanned aerial systems for precision agriculture: A review. Precis Agric 2012; 13(6): 693-712. [http://dx.doi.org/10.1007/s11119-012-9274-5]",
  weiss: "Weiss M, Jacob F, Duveiller G. Remote sensing for agricultural applications: A meta-review. Remote Sens Environ 2020; 236: 111402. [http://dx.doi.org/10.1016/j.rse.2019.111402]",
  gorelick: "Gorelick N, Hancher M, Dixon M, Ilyushchenko S, Thau D, Moore R. Google Earth Engine: Planetary-scale geospatial analysis for everyone. Remote Sens Environ 2017; 202: 18-27. [http://dx.doi.org/10.1016/j.rse.2017.06.031]",
  sishodia: "Sishodia RP, Ray RL, Singh SK. Applications of remote sensing in precision agriculture: A review. Remote Sens (Basel) 2020; 12(19): 3136. [http://dx.doi.org/10.3390/rs12193136]",
  mohanty: "Mohanty SP, Hughes DP, Salathé M. Using deep learning for image-based plant disease detection. Front Plant Sci 2016; 7: 1419. [http://dx.doi.org/10.3389/fpls.2016.01419]",
  hughes: "Hughes DP, Salathé M. An open access repository of images on plant health to enable the development of mobile disease diagnostics. arXiv 2015; 1511.08060.",
  ferentinos: "Ferentinos KP. Deep learning models for plant disease detection and diagnosis. Comput Electron Agric 2018; 145: 311-8. [http://dx.doi.org/10.1016/j.compag.2018.01.009]",
  barbedo: "Barbedo JGA. Factors influencing the use of deep learning for plant disease recognition. Biosyst Eng 2018; 172: 84-91. [http://dx.doi.org/10.1016/j.biosystemseng.2018.05.013]",
  mahlein: "Mahlein AK. Plant disease detection by imaging sensors – Parallels and specific demands for precision agriculture and plant phenotyping. Plant Dis 2016; 100(2): 241-51. [http://dx.doi.org/10.1094/PDIS-03-15-0340-FE]",
  goap: "Goap A, Sharma D, Shukla AK, Rama Krishna C. An IoT based smart irrigation management system using machine learning and open source technologies. Comput Electron Agric 2018; 155: 41-9. [http://dx.doi.org/10.1016/j.compag.2018.09.040]",
  pantazi: "Pantazi XE, Moshou D, Alexandridis T, Whetton RL, Mouazen AM. Wheat yield prediction using machine learning and advanced sensing techniques. Comput Electron Agric 2016; 121: 57-65. [http://dx.doi.org/10.1016/j.compag.2015.11.018]",
  khaki: "Khaki S, Wang L. Crop yield prediction using deep neural networks. Front Plant Sci 2019; 10: 621. [http://dx.doi.org/10.3389/fpls.2019.00621]",
  vanklomp: "van Klompenburg T, Kassahun A, Catal C. Crop yield prediction using machine learning: A systematic literature review. Comput Electron Agric 2020; 177: 105709. [http://dx.doi.org/10.1016/j.compag.2020.105709]",
  bechar: "Bechar A, Vigneault C. Agricultural robots for field operations: Concepts and components. Biosyst Eng 2016; 149: 94-111. [http://dx.doi.org/10.1016/j.biosystemseng.2016.06.014]",
  shamshiri: "Shamshiri RR, Weltzien C, Hameed IA, et al. Research and development in agricultural robotics: A perspective of digital farming. Int J Agric Biol Eng 2018; 11(4): 1-14. [http://dx.doi.org/10.25165/j.ijabe.20181104.4278]",
  lowenberg: "Lowenberg-DeBoer J, Huang IY, Grigoriadis V, Blackmore S. Economics of robots and automation in field crop production. Precis Agric 2020; 21(2): 278-99. [http://dx.doi.org/10.1007/s11119-019-09667-5]",
  sharma: "Sharma R, Kamble SS, Gunasekaran A, Kumar V, Kumar A. A systematic literature review on machine learning applications for sustainable agriculture supply chain performance. Comput Oper Res 2020; 119: 104926. [http://dx.doi.org/10.1016/j.cor.2020.104926]",
  arrieta: "Barredo Arrieta A, Díaz-Rodríguez N, Del Ser J, et al. Explainable Artificial Intelligence (XAI): Concepts, taxonomies, opportunities and challenges toward responsible AI. Inf Fusion 2020; 58: 82-115. [http://dx.doi.org/10.1016/j.inffus.2019.12.012]",
  lundberg: "Lundberg SM, Lee SI. A unified approach to interpreting model predictions. Adv Neural Inf Process Syst 2017; 30: 4765-74.",
  ryan: "Ryan M. Agricultural big data analytics and the ethics of power. J Agric Environ Ethics 2020; 33(1): 49-69. [http://dx.doi.org/10.1007/s10806-019-09812-0]",
  mcmahan: "McMahan B, Moore E, Ramage D, Hampson S, Agüera y Arcas B. Communication-efficient learning of deep networks from decentralized data. Proceedings of the 20th International Conference on Artificial Intelligence and Statistics (AISTATS). Fort Lauderdale, FL, USA. 2017; 54: pp. 1273-82.",
};

const order = [];
function num(key) {
  if (!REFS[key]) throw new Error("Unknown reference key: " + key);
  let i = order.indexOf(key);
  if (i < 0) { order.push(key); i = order.length - 1; }
  return i + 1;
}
function cite(keys) {
  const ns = [...new Set(keys.map(num))].sort((a, b) => a - b);
  const parts = [];
  for (let i = 0; i < ns.length; ) {
    let j = i;
    while (j + 1 < ns.length && ns[j + 1] === ns[j] + 1) j++;
    if (j - i >= 2) parts.push(`${ns[i]}-${ns[j]}`);
    else for (let k = i; k <= j; k++) parts.push(String(ns[k]));
    i = j + 1;
  }
  return "[" + parts.join(", ") + "]";
}

// Inline markup: {key1,key2} -> citation; *italic*; **bold**
function runs(text, base = {}) {
  text = text.replace(/\{([a-zA-Z0-9_,\s]+)\}/g, (_, k) => cite(k.split(",").map((s) => s.trim())));
  const out = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*)/g;
  let last = 0, m;
  while ((m = re.exec(text))) {
    if (m.index > last) out.push(new TextRun({ text: text.slice(last, m.index), font: FONT, ...base }));
    const t = m[0];
    if (t.startsWith("**")) out.push(new TextRun({ text: t.slice(2, -2), bold: true, font: FONT, ...base }));
    else out.push(new TextRun({ text: t.slice(1, -1), italics: true, font: FONT, ...base }));
    last = m.index + t.length;
  }
  if (last < text.length) out.push(new TextRun({ text: text.slice(last), font: FONT, ...base }));
  return out;
}

const body = [];
const LINE = 360; // 1.5 line spacing
const P = (t) => body.push(new Paragraph({ children: runs(t), alignment: AlignmentType.JUSTIFIED, spacing: { line: LINE, after: 160 }, indent: { firstLine: 0 } }));
const H1 = (t) => body.push(new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun({ text: t.toUpperCase(), bold: true, font: FONT, size: 24 })], spacing: { before: 320, after: 160, line: LINE } }));
const H2 = (t) => body.push(new Paragraph({ heading: HeadingLevel.HEADING_2, children: [new TextRun({ text: t, bold: true, font: FONT, size: 24 })], spacing: { before: 240, after: 120, line: LINE } }));
const H3 = (t) => body.push(new Paragraph({ heading: HeadingLevel.HEADING_3, children: [new TextRun({ text: t, bold: true, italics: true, font: FONT, size: 24 })], spacing: { before: 200, after: 100, line: LINE } }));
const BUL = (t) => body.push(new Paragraph({ numbering: { reference: "bullets", level: 0 }, children: runs(t), alignment: AlignmentType.JUSTIFIED, spacing: { line: LINE, after: 80 } }));

function FIGURE(file, n, caption, wpx, hpx) {
  body.push(new Paragraph({
    alignment: AlignmentType.CENTER, spacing: { before: 200, after: 80 },
    children: [new ImageRun({ type: "png", data: fs.readFileSync(path.join(FIG, file)), transformation: { width: wpx, height: hpx },
      altText: { title: `Fig. (${n})`, description: caption, name: file } })],
  }));
  body.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 240 },
    children: [new TextRun({ text: `Fig. (${n}). `, bold: true, font: FONT, size: 22 }), ...runs(caption, { size: 22 })] }));
}

const TABLE_W = 9026; // A4 width minus 1" margins (DXA)
const border = { style: BorderStyle.SINGLE, size: 4, color: "808080" };
const borders = { top: border, bottom: border, left: border, right: border };
function TABLE(n, caption, header, rows, widths) {
  body.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 240, after: 120 }, keepNext: true,
    children: [new TextRun({ text: `Table ${n}. `, bold: true, font: FONT, size: 22 }), ...runs(caption, { size: 22 })] }));
  const cell = (t, w, hdr) => new TableCell({
    borders, width: { size: w, type: WidthType.DXA },
    shading: hdr ? { fill: "D9EAD3", type: ShadingType.CLEAR, color: "auto" } : undefined,
    margins: { top: 60, bottom: 60, left: 100, right: 100 },
    children: [new Paragraph({ children: runs(t, { size: 20, bold: hdr || undefined }) })],
  });
  body.push(new Table({
    width: { size: TABLE_W, type: WidthType.DXA }, columnWidths: widths,
    rows: [new TableRow({ tableHeader: true, children: header.map((h, i) => cell(h, widths[i], true)) }),
      ...rows.map((r) => new TableRow({ children: r.map((c, i) => cell(c, widths[i], false)) }))],
  }));
  body.push(new Paragraph({ spacing: { after: 200 }, children: [] }));
}

// ---------------------------------------------------------------------------
// CHAPTER CONTENT
// ---------------------------------------------------------------------------
H1("Introduction");
P("Agriculture is the oldest organised human enterprise and remains the foundation of food security, rural livelihoods and national economies. It is now under pressure from several directions at once. Global population is expected to approach ten billion by the middle of the century, dietary patterns are shifting towards more resource-intensive foods, and the natural resources on which farming depends (fertile soil, fresh water, biodiversity and a stable climate) are under growing stress {godfray,un2022}. The Food and Agriculture Organization of the United Nations (FAO) has estimated that overall food production would need to rise by roughly 70 per cent between the late 2000s and 2050 to meet projected demand {fao2009}. This has to happen on a largely fixed land base and with a smaller environmental footprint {tilman}.");
P("Mechanisation, improved crop varieties, irrigation and agrochemicals greatly increased productivity during the twentieth century. However, they relied on uniform management of fields that are in fact highly variable in space and time. Large gains from blanket application of inputs are now harder to find, and their environmental costs are well documented. The central question for modern agriculture has therefore changed from *how to produce more* to *how to produce more with less, and more intelligently*. Making better decisions, at the right place, at the right time and at the right scale, has become the key lever of agricultural progress.");
P("Artificial Intelligence (AI) is well suited to this task. AI is the broad discipline of building computational systems that perceive their environment, learn from data, reason under uncertainty and act towards goals. Several developments have converged to make AI practical in agriculture: low-cost sensors, ubiquitous connectivity, satellite and drone imagery, cloud computing and, above all, machine learning and deep learning. Together they allow a farm to be observed continuously and at fine resolution. Its condition can then be interpreted automatically, and the resulting insight can be turned into specific actions {liakos,kamilaris,jha}.");
P("This introductory chapter lays the conceptual groundwork for the rest of the book. It first examines why agriculture needs AI, by reviewing the global challenges facing food systems. It then traces the evolution of agricultural technology from mechanisation to the data-driven paradigms known as Agriculture 4.0 and 5.0. Next, it introduces the fundamental concepts of AI and machine learning in terms accessible to readers from both the computational and the agricultural sciences. A reference architecture for AI-enabled agricultural systems follows, along with a survey of the main application domains and their benefits. The chapter then discusses the technical, economic and ethical challenges that currently limit adoption. It closes with emerging directions and a roadmap of the chapters that follow.");

H1("The Global Agricultural Challenge and the Case for AI");
H2("Food Security and Sustainable Intensification");
P("Food security has four dimensions: availability, access, utilisation and stability. All four are sensitive to how efficiently and resiliently food is produced. Godfray *et al.* argued that meeting future demand will require action on several fronts simultaneously: closing yield gaps, increasing production limits, reducing waste, changing diets and expanding aquaculture {godfray}. Tilman *et al.* showed that the environmental consequences of meeting demand depend strongly on the pathway chosen. Raising yields on existing farmland, especially in low-yielding regions, could meet demand with far less land clearing and greenhouse-gas emission than expanding cropland {tilman}. This strategy of *sustainable intensification* is inherently information-intensive. It requires knowing where, when and how much to intervene.");
P("The structure of world agriculture adds a further dimension. Of the more than 570 million farms worldwide, the large majority are small, family-operated holdings, most of them in low- and middle-income countries {lowder}. Any technological transformation that aims to have a global impact must therefore work not only for large, capital-intensive farms but also for smallholders with limited resources, connectivity and formal training. This theme recurs throughout the chapter.");

H2("Climate, Resource and Labour Constraints");
P("Climate change increases the variability of temperature and rainfall and the frequency of extreme events such as droughts, heatwaves and floods. It also shifts the geographic range of pests and pathogens. Under these conditions, historical experience and fixed calendars become less reliable guides to farm management. At the same time, agriculture is the largest consumer of fresh water, and excessive or poorly timed fertiliser and pesticide use pollutes soil and water. Many regions also face shortages of agricultural labour as rural populations age and migrate to cities. These pressures create a demand for technologies that can sense changing conditions in real time, anticipate risks and automate labour-intensive tasks.");

H2("Limits of Conventional Decision-Making");
P("Traditional farm management relies on the experience of the farmer, periodic manual scouting and generalised recommendations from extension services. This knowledge is valuable but has inherent limits. Human scouting can sample only a small fraction of a field. Visual symptoms of stress often appear only after yield damage has begun. Recommendations calibrated for an average field ignore within-field variability in soil, topography and microclimate {gebbers}. Agricultural systems are also complex and non-linear: yield emerges from interactions among genotype, environment and management that are difficult to capture with simple rules. AI and data-driven methods are attractive because they can integrate large, heterogeneous data streams and learn these complex relationships directly from observations {wolfert,liakos}.");

H1("Evolution of Technology in Agriculture");
P("The adoption of AI in agriculture is best understood as the latest stage in a long sequence of technological revolutions. Each stage has increased the quantity of information available to decision-makers and the degree to which actions can be automated. These stages are commonly described as Agriculture 1.0 to 5.0, as summarised in Fig. (**1**).");
FIGURE("Fig1_evolution.png", 1, "Evolution of agricultural technology from Agriculture 1.0 to Agriculture 5.0, showing the progressive increase in data intensity, connectivity and decision autonomy.", 620, 223);

H2("From Mechanisation to Precision Agriculture");
P("Agriculture 1.0 describes traditional farming based on human and animal labour. Agriculture 2.0 began with mechanisation, the internal combustion tractor and, in the mid-twentieth century, the Green Revolution. The Green Revolution combined high-yielding varieties with synthetic fertilisers, pesticides and irrigation. These advances greatly increased output but encouraged uniform treatment of large areas.");
P("Agriculture 3.0, commonly called precision agriculture, emerged in the 1980s and 1990s. It was enabled by satellite positioning (GNSS), yield monitors, geographic information systems and variable-rate application technology {zhangN,gebbers}. Its core principle is to manage spatial and temporal variability by applying the right input at the right rate, place and time. Remote sensing from satellites and aircraft provided field-scale maps of crop vigour, and vegetation indices became standard diagnostic tools {mulla}. However, the analysis of these data remained largely manual or based on simple statistical models.");

H2("Digital Agriculture and Agriculture 4.0");
P("Agriculture 4.0, also called digital agriculture or smart farming, uses the Internet of Things (IoT), cloud computing and big-data analytics to connect sensors, machines and people across the farm and the wider value chain {wolfert,klerkx}. Data are no longer collected occasionally for a single purpose. Instead they stream continuously from soil probes, weather stations, cameras, drones, machinery and market platforms. The resulting volume, velocity and variety of data exceed what can be interpreted manually. Machine learning therefore becomes an essential tool rather than an optional one {liakos,basso}.");

H2("Towards Agriculture 5.0");
P("Agriculture 5.0 is the emerging stage in which AI moves from providing insight to taking autonomous or semi-autonomous action. Examples include robots that weed or harvest selectively, irrigation systems that schedule themselves and decision-support platforms that continuously learn from outcomes. Digital twins, which are virtual replicas of fields, animals or entire farms that are synchronised with real-time data, allow scenarios to be simulated before they are applied {verdouw}. Consistent with the broader Industry 5.0 vision, Agriculture 5.0 also emphasises human-centricity, resilience and sustainability. It aims to augment rather than replace the farmer and requires innovation to be developed responsibly with stakeholders {rose}.");

H1("Fundamental Concepts of Artificial Intelligence");
H2("Defining Artificial Intelligence");
P("The term *artificial intelligence* was introduced in the 1955 proposal for the Dartmouth summer research project. That proposal conjectured that every aspect of learning or intelligence could in principle be described precisely enough for a machine to simulate it {mccarthy}. In contemporary textbooks AI is framed as the study of *rational agents*: systems that perceive their environment through sensors and act on it through actuators so as to achieve the best expected outcome {russell}. This agent view maps naturally onto agriculture. The environment is the field, orchard, greenhouse or herd. Sensors include cameras, probes and satellites. Actuators range from sprayers and valves to robots and, very often, human farmers acting on recommendations.");
P("AI encompasses several families of techniques. These are summarised in Fig. (**2**) and Table **1** and described briefly below. Chapter 2 treats their scientific foundations in depth.");

H2("Knowledge-Based Systems and Fuzzy Logic");
P("The first generation of agricultural AI consisted of expert systems. These encode the knowledge of human specialists as explicit rules, for example *IF leaf lesions are concentric AND humidity is high THEN suspect early blight*. Such systems were developed from the 1980s for crop management, pest diagnosis and irrigation advice {mckinion}. Fuzzy logic, which handles imprecise linguistic concepts such as *dry*, *warm* or *severe* through degrees of membership rather than binary truth values {zadeh}, remains widely used in irrigation and greenhouse climate controllers. Knowledge-based approaches are transparent and need little data. However, they are laborious to build and maintain, and they cannot easily discover patterns that experts have not already articulated.");

H2("Machine Learning");
P("Machine Learning (ML) reverses the knowledge-engineering approach. Instead of hand-crafting rules, an algorithm learns a mapping from data. ML is usually divided into three paradigms:");
BUL("**Supervised learning** learns from labelled examples, such as images tagged with disease names or historical records pairing weather and soil data with yields. It covers *classification*, which predicts a category, and *regression*, which predicts a continuous quantity. Widely used algorithms include support vector machines {cortes}, decision trees and random forests {breiman}, gradient boosting and artificial neural networks.");
BUL("**Unsupervised learning** finds structure in unlabelled data. Examples include clustering fields into management zones, detecting anomalous sensor readings and reducing the dimensionality of hyperspectral images.");
BUL("**Reinforcement learning** trains an agent to choose actions that maximise cumulative reward through interaction with an environment {sutton}. It is well suited to sequential decision problems such as irrigation scheduling, greenhouse climate control and robot navigation, typically trained first in simulation.");
P("Liakos *et al.* reviewed ML applications across crop, livestock, water and soil management. They found that supervised methods dominated the literature, with artificial neural networks and support vector machines among the most frequently used models {liakos}.");

H2("Deep Learning");
P("Deep Learning (DL) uses neural networks with many layers that learn hierarchical representations directly from raw data. This removes the need for manually engineered features {lecun}. The breakthrough of convolutional neural networks (CNNs) in large-scale image classification in 2012 {krizhevsky}, followed by deeper architectures such as residual networks {he2016}, transformed computer vision. Real-time object detectors such as YOLO {redmon} made it feasible to locate individual plants, weeds, fruits or animals in video streams. For sequential data such as weather series and sensor logs, recurrent networks and long short-term memory (LSTM) units {hochreiter} and, more recently, attention-based Transformers {vaswani} are widely used. Kamilaris and Prenafeta-Boldú surveyed 40 early applications of DL in agriculture and found that it generally outperformed conventional image-processing and ML approaches, while noting its dependence on large annotated datasets {kamilaris}.");

H2("Foundation and Generative Models");
P("The most recent development is *foundation models*: very large models pre-trained on broad data that can be adapted to many downstream tasks with little additional labelled data. Examples include large language models capable of conversational question answering and vision models such as the Segment Anything Model, which produces object masks for previously unseen images {kirillov}. In agriculture, foundation models are being explored for multilingual farmer advisory chatbots, zero-shot or few-shot crop and pest recognition, and synthetic data generation to augment scarce training sets {li2024}. They promise to lower the data barrier that has limited specialised models. However, they raise new questions of accuracy, domain grounding, computational cost and accountability.");
FIGURE("Fig2_taxonomy.png", 2, "Taxonomy of AI techniques used in agriculture, showing representative methods for each family and typical agricultural tasks to which they are applied.", 620, 335);
TABLE(1, "Major AI paradigms, their principles and representative agricultural applications.",
  ["AI paradigm", "Underlying principle", "Strengths", "Limitations", "Representative applications"],
  [
    ["Expert systems and fuzzy logic", "Explicit rules and degrees of membership elicited from domain experts", "Transparent; little data needed", "Costly to build and update; limited adaptability", "Pest diagnosis advisory, greenhouse and irrigation controllers"],
    ["Classical machine learning (SVM, RF, gradient boosting)", "Statistical learning of mappings from engineered features", "Effective on small and tabular datasets; relatively interpretable", "Depends on feature engineering; limited on raw images", "Yield prediction, soil property estimation, management-zone delineation"],
    ["Deep learning (CNN, LSTM, Transformer)", "Hierarchical representation learning from raw data", "State-of-the-art accuracy on images and sequences", "Data-hungry; computationally intensive; opaque", "Disease and weed detection, fruit counting, time-series forecasting"],
    ["Reinforcement learning", "Learning action policies from rewards through interaction", "Handles sequential decisions under uncertainty", "Needs simulators; safety of exploration", "Irrigation scheduling, climate control, robot path planning"],
    ["Foundation and generative models", "Large-scale pre-training followed by adaptation", "Few-shot generalisation; natural-language interaction", "Hallucination risk; high resource cost; domain grounding", "Farmer advisory chatbots, synthetic training data, zero-shot segmentation"],
  ],
  [1700, 1900, 1750, 1750, 1926]);

H1("Architecture of an AI-Enabled Agricultural System");
P("An AI application in agriculture is never an algorithm in isolation. It is embedded in a socio-technical system that spans physical sensing, communication, data management, analytics and action. Fig. (**3**) presents a five-layer reference architecture synthesised from the smart-farming literature {wolfert,tzounis,elijah}. Cross-cutting concerns of security, privacy, explainability and ethics apply to every layer.");
FIGURE("Fig3_architecture.png", 3, "Layered reference architecture of an AI-enabled agricultural system, with cross-cutting concerns of security, privacy, explainability and ethics.", 520, 403);

H2("Sensing and Data Acquisition");
P("The perception layer captures the state of the crop, soil, animals and environment. In-field sensors measure soil moisture, temperature, electrical conductivity, nutrient levels and micro-weather {tzounis}. Cameras operating in visible, multispectral, hyperspectral and thermal bands reveal plant vigour, water stress and disease symptoms. Unmanned aerial vehicles (UAVs) provide very high-resolution imagery on demand {zhangK}, while Earth-observation satellites offer frequent, wide-area coverage over many seasons {weiss,sishodia}. Modern machinery also generates telemetry, yield maps and application records. Farmer observations and market data complete the picture. Chapter 3 and Chapter 10 examine these acquisition technologies in detail.");

H2("Connectivity and Edge Computing");
P("Data must be moved from often remote fields to where they can be processed. Low-power wide-area networks such as LoRaWAN and NB-IoT suit battery-powered sensors that send small packets over long distances. Cellular networks, including 5G, support high-bandwidth applications such as video. Because rural connectivity is frequently intermittent, *edge computing* is increasingly important. In this approach, inference runs directly on gateways, drones, tractors or microcontrollers (TinyML), which reduces latency and bandwidth and allows operation offline {elijah}. Hardware design for low-power edge inference is an active area where advances in embedded systems and VLSI directly benefit agricultural AI.");

H2("Data Management");
P("Raw agricultural data are heterogeneous, noisy and incomplete. This layer performs cleaning, calibration, geo-referencing and temporal alignment. It also fuses data from multiple sources, annotates training data, and stores everything with appropriate metadata in cloud or on-premise repositories. Cloud geospatial platforms such as Google Earth Engine have made petabyte-scale satellite archives accessible for analysis without local infrastructure {gorelick}. Data governance, which defines who owns, may access and may reuse farm data, is a crucial and often neglected element of this layer {wolfert}.");

H2("AI and Analytics");
P("The analytics layer hosts the models introduced earlier in this chapter. Examples include classifiers that recognise diseases, detectors that count fruit, regressors that forecast yield, optimisers that compute input prescriptions and knowledge graphs that link agronomic concepts. Models may be trained centrally and deployed to the edge, or trained collaboratively across farms without pooling raw data through *federated learning* {mcmahan}. This approach helps preserve the privacy of commercially sensitive farm records.");

H2("Decision Support and Actuation");
P("Finally, insights must lead to action. This can occur through dashboards and mobile advisories that inform human decisions, or through prescription maps executed by variable-rate machinery. It can also take the form of fully autonomous actuation, for example a robot that removes a weed or a valve that opens an irrigation zone. The choice of automation level depends on the reliability of the model, the cost of errors and the preferences of the farmer.");

H2("The Closed-Loop Learning Cycle");
P("These layers operate as a continuous cycle rather than a one-way pipeline (Fig. **4**). Each action changes the state of the field, and the consequences are observed in the next round of sensing. These outcomes provide feedback that can be used to evaluate and retrain models. Over successive seasons, this loop lets an AI-enabled farm improve its decisions by learning from its own experience. A digital twin formalises the loop by maintaining a continuously updated virtual counterpart of the physical system {verdouw}.");
FIGURE("Fig4_loop.png", 4, "The closed-loop observe–decide–act–learn cycle that underlies AI-enabled agricultural management.", 420, 420);

H1("Major Application Domains");
P("AI is applied across the entire agricultural value chain, from pre-season planning to the consumer's table. Table **2** maps representative applications to value-chain stages. The following subsections introduce the principal domains, each of which is developed in a dedicated chapter of this book.");
TABLE(2, "Representative AI applications across the agricultural value chain.",
  ["Value-chain stage", "Decision or task", "Typical data sources", "Typical AI techniques"],
  [
    ["Pre-production / planning", "Crop and variety selection, land suitability, soil mapping", "Soil surveys, historical weather, satellite archives", "Random forest, gradient boosting, geostatistics with ML"],
    ["Crop establishment and growth", "Crop-health monitoring, nutrient and irrigation management", "IoT soil sensors, UAV and satellite imagery, weather stations", "CNN, vegetation-index regression, reinforcement learning, fuzzy control"],
    ["Crop protection", "Detection of diseases, pests and weeds; targeted spraying", "RGB, multispectral and hyperspectral images; trap images", "CNN classifiers, object detectors (e.g., YOLO), segmentation networks"],
    ["Harvest", "Yield estimation, maturity assessment, robotic harvesting", "Field images, yield monitors, machine telemetry", "Object detection and counting, LSTM, robot perception and control"],
    ["Livestock and aquaculture", "Health and behaviour monitoring, feeding optimisation", "Wearable sensors, video, water-quality sensors", "Time-series classification, computer vision, anomaly detection"],
    ["Post-harvest and supply chain", "Grading, storage monitoring, demand and price forecasting, traceability", "Machine-vision lines, storage sensors, market data", "CNN grading, time-series forecasting, optimisation"],
  ],
  [1800, 2500, 2400, 2326]);

H2("Crop and Soil Monitoring");
P("Continuous monitoring of crop and soil condition underpins nearly every other application. ML models translate spectral reflectance from satellites and drones into estimates of biophysical variables such as leaf area index, chlorophyll content, biomass and canopy water status {weiss,sishodia}. They also estimate soil properties such as organic carbon and texture. Time series of such estimates reveal crop phenology and detect stress before it is visible to the eye, so that intervention can be timely and site-specific {mulla}.");

H2("Plant Disease, Pest and Weed Detection");
P("Image-based diagnosis is the most widely studied application of deep learning in agriculture. Mohanty *et al.* trained CNNs on the public PlantVillage dataset of 54,306 leaf images covering 14 crop species and 26 diseases {hughes}. They reported accuracies above 99 per cent on held-out images, but performance dropped markedly on images captured under conditions different from the training set {mohanty}. Ferentinos subsequently reported similarly high accuracies on a larger dataset of 87,848 images spanning 25 plant species {ferentinos}. These studies demonstrate both the promise of the approach and its central weakness. Models trained on clean, homogeneous images often generalise poorly to real field conditions with complex backgrounds, variable lighting and multiple co-occurring symptoms {barbedo}. Hyperspectral and thermal sensing can detect infection before visible symptoms appear {mahlein}. Object detectors enable site-specific weeding and precision spraying. These topics are the subject of Chapters 4, 5 and 7.");

H2("Precision Irrigation and Nutrient Management");
P("AI-driven irrigation combines soil-moisture measurements, weather forecasts and crop models to predict water requirements and schedule irrigation automatically. Goap *et al.*, for example, combined IoT sensing with ML-based prediction of soil moisture and open-source technologies to build a low-cost smart irrigation system {goap}. Similar approaches support variable-rate fertiliser application that matches nutrient supply to spatially varying crop demand, reducing both cost and environmental losses {gebbers}. Chapters 6 and 7 develop these themes.");

H2("Yield Prediction and Agricultural Decision-Making");
P("Accurate, early yield forecasts support decisions by farmers, insurers, traders and governments. ML models that fuse soil, weather, management and remote-sensing data have been shown to predict yield with useful accuracy at field and regional scales {pantazi,khaki}. A systematic review by van Klompenburg *et al.* found that temperature, rainfall and soil type were the most common input features, and that artificial neural networks were the most frequently applied algorithms, with CNNs and LSTMs being the most common deep-learning architectures {vanklomp}. Beyond yield, predictive analytics support planting-date selection, risk assessment and market planning (Chapter 9).");

H2("Agricultural Robotics and Automation");
P("Agricultural robots combine AI-based perception with navigation and manipulation to perform tasks such as weeding, spraying, pruning, phenotyping and harvesting {bechar,shamshiri}. Unstructured, changing outdoor environments and the delicacy of biological objects make agricultural robotics more challenging than factory automation. Recent advances in computer vision, sensor fusion and motion planning are steadily expanding what is feasible. Economic analyses suggest that small autonomous machines could change the economics of field operations and make precision management viable on smaller fields {lowenberg}. Chapter 8 addresses autonomous agricultural robots and intelligent field systems.");

H2("Livestock and Aquaculture");
P("In animal production, AI analyses data from wearable accelerometers, cameras, microphones and milking systems. These analyses detect lameness, oestrus and disease early, monitor feeding behaviour and assess welfare {liakos}. In aquaculture, ML models predict water-quality parameters and optimise feeding to reduce waste. Although this book concentrates on crop systems, the same principles of sensing, learning and closed-loop control apply.");

H2("Supply Chain, Markets and Post-Harvest");
P("Beyond the farm gate, machine-vision systems grade and sort produce, and sensors with predictive models monitor storage conditions to reduce post-harvest losses. Forecasting models anticipate demand and prices, while optimisation algorithms improve logistics. Sharma *et al.* reviewed ML applications across the agricultural supply chain and highlighted their potential to improve sustainability and performance from production planning to distribution {sharma}.");

H1("Benefits and Impacts");
P("When appropriately designed and deployed, AI can deliver benefits across the economic, environmental and social dimensions of sustainability {liakos,jha,basso}:");
BUL("**Productivity:** earlier detection of stress and more timely, better-targeted interventions help close yield gaps.");
BUL("**Resource efficiency:** site-specific application of water, fertiliser and pesticides reduces input costs and losses to the environment.");
BUL("**Risk reduction and resilience:** forecasts and early warnings of pests, diseases and weather extremes help farmers and institutions prepare and respond.");
BUL("**Labour relief:** automation of repetitive and physically demanding tasks addresses labour shortages and improves working conditions.");
BUL("**Knowledge access:** digital advisory services can bring expert-level diagnosis and recommendations to farmers who lack access to extension services.");
BUL("**Traceability and transparency:** data-driven records support quality assurance, certification and consumer trust.");
P("These benefits are not automatic. They depend on the accuracy and robustness of the models, the quality of the underlying data, the appropriateness of the technology for the local context and the extent to which farmers trust and adopt it. The following section examines these conditions.");

H1("Challenges and Limitations");
P("Despite rapid research progress, the use of AI on farms, particularly smallholder farms, remains limited. Table **3** summarises the main barriers and possible mitigation strategies, which are discussed below.");
TABLE(3, "Key challenges for AI adoption in agriculture and possible mitigation strategies.",
  ["Challenge", "Description", "Possible mitigation strategies"],
  [
    ["Data scarcity and quality", "Few large, well-annotated, representative datasets; labels costly and requiring expert knowledge", "Open data initiatives, transfer and self-supervised learning, data augmentation, citizen-science labelling"],
    ["Poor generalisation", "Models trained in one region, season or sensor degrade elsewhere (domain shift)", "Multi-site training, domain adaptation, field-realistic benchmarks, continual learning"],
    ["Infrastructure and cost", "Limited rural connectivity, power and capital; high cost of sensors and machinery", "Edge/offline inference, low-cost open-source hardware, service-based business models"],
    ["Digital literacy and the digital divide", "Uneven skills and access risk widening inequality between farmers", "Participatory design, local-language interfaces, extension and training programmes"],
    ["Explainability and trust", "Opaque models hinder acceptance and accountability", "Explainable AI methods, uncertainty quantification, human-in-the-loop design"],
    ["Data ownership, privacy and ethics", "Unclear rights over farm data; power concentration; algorithmic bias", "Data-sharing codes of conduct, federated learning, responsible-innovation frameworks"],
  ],
  [2000, 3400, 3626]);

H2("Data Availability and Quality");
P("Supervised AI requires large volumes of labelled data that are representative of the conditions in which the model will be used. Agricultural data are expensive to collect, seasonal (often only one growing cycle per year) and highly diverse across crops, varieties, regions and growth stages. Labelling frequently requires agronomic or pathological expertise. Many published datasets were collected under controlled conditions and do not capture the complexity of real fields {kamilaris,barbedo}.");

H2("Generalisation and Robustness");
P("A model that performs well on its test set may fail when deployed in a different region, season or with a different camera. This phenomenon is known as *domain shift*. The marked drop in accuracy observed when plant-disease models trained on laboratory images were applied to field images is a well-known example {mohanty,barbedo}. Robust evaluation on independent, field-realistic data from multiple sites and seasons is therefore essential before deployment claims can be made.");

H2("Infrastructure, Cost and the Digital Divide");
P("Many rural areas lack reliable internet connectivity and electricity. Precision equipment, sensors and robots require substantial investment that may not be justified on small or fragmented holdings. Unless these barriers are addressed, AI risks benefiting mainly large, well-capitalised farms and widening existing inequalities {klerkx,lowder}. Low-cost, open-source hardware and software, edge inference and shared-service models such as custom-hiring centres or drone-as-a-service are promising ways to widen access.");

H2("Explainability and Trust");
P("Farmers are understandably reluctant to act on recommendations they cannot understand, especially when their livelihood is at stake. Deep neural networks are often regarded as *black boxes*. Explainable AI (XAI) provides methods to make model behaviour interpretable, such as feature-attribution techniques like SHAP {lundberg} and saliency maps that highlight the image regions driving a diagnosis {arrieta}. Quantifying and communicating the uncertainty of predictions is equally important. Chapter 11 is devoted to explainable, trustworthy and responsible AI in agriculture.");

H2("Data Ownership, Privacy and Ethics");
P("Farm data have considerable commercial value, and questions about who owns, controls and profits from them are not yet settled. Concentration of data and analytical capability among a small number of technology and agribusiness firms may shift power away from farmers {ryan,wolfert}. Other ethical concerns include algorithmic bias against under-represented crops or regions, the effects of automation on employment and accountability when automated decisions cause harm. Responsible innovation frameworks call for anticipating such impacts and involving farmers and other stakeholders throughout the design process {rose}.");

H1("Future Directions");
P("Several converging trends are likely to shape the next decade of AI in agriculture. Chapter 12 explores them in greater depth:");
BUL("**Foundation models and generative AI** adapted to agriculture, enabling multilingual conversational advisory services and few-shot recognition of new crops and pests {li2024}.");
BUL("**Digital twins** of fields, farms and supply chains that integrate process-based crop models with data-driven learning for scenario analysis and optimisation {verdouw}.");
BUL("**Edge AI and energy-efficient hardware**, including low-power AI accelerators and TinyML on microcontrollers. These will bring intelligence to battery-powered sensors and small robots operating offline.");
BUL("**Federated and privacy-preserving learning** that allows farms and organisations to build shared models without surrendering their raw data {mcmahan}.");
BUL("**Swarm and collaborative robotics**, in which fleets of small, lightweight robots and drones cooperate to perform field operations with less soil compaction {lowenberg}.");
BUL("**Hybrid physics-informed AI** that combines agronomic knowledge and mechanistic models with ML, improving data efficiency, interpretability and robustness under climate change.");
BUL("**Responsible and inclusive AI**, with emphasis on fairness, accessibility for smallholders, local-language interfaces and co-design with farming communities {rose,klerkx}.");
P("Realising these opportunities will require sustained interdisciplinary collaboration. Agronomists, plant pathologists and soil scientists need to work alongside computer scientists, electronics and communication engineers, roboticists, economists and social scientists, and above all with farmers themselves.");

H1("Organisation of this Book");
P("This book is organised to lead the reader from foundations through enabling technologies and specific applications to cross-cutting issues and future prospects. Table **4** summarises the scope of each chapter.");
TABLE(4, "Organisation of the book and scope of each chapter.",
  ["Chapter", "Title", "Scope"],
  [
    ["1", "Introduction to Artificial Intelligence in Agriculture", "Motivation, technological evolution, AI concepts, system architecture, applications, challenges and outlook"],
    ["2", "Scientific Foundations of AI and Machine Learning for Agriculture", "Mathematical and algorithmic basis of ML and DL models"],
    ["3", "Agricultural Data Acquisition and IoT-Based Sensing Systems", "Sensors, networks and data pipelines"],
    ["4", "Computer Vision and Deep Learning for Crop Monitoring", "Image analysis for crop growth and health"],
    ["5", "AI for Plant Disease and Pest Detection", "Diagnosis, early warning and surveillance"],
    ["6", "AI-Driven Precision Irrigation and Water Management", "Water-demand prediction and smart irrigation control"],
    ["7", "AI for Precision Fertilizer and Pesticide Application", "Variable-rate and targeted input application"],
    ["8", "Autonomous Agricultural Robots and Intelligent Field Systems", "Perception, navigation and manipulation in the field"],
    ["9", "Predictive Analytics for Crop Yield and Agricultural Decision-Making", "Forecasting and decision support"],
    ["10", "Remote Sensing, Drones, and Satellite-Based AI Applications", "UAV and Earth-observation analytics"],
    ["11", "Explainable, Trustworthy, and Responsible AI in Agriculture", "Interpretability, ethics, governance and trust"],
    ["12", "Future Directions and Emerging Trends in AI-Enabled Agriculture", "Emerging technologies and research agenda"],
  ],
  [1000, 4000, 4026]);

H1("Conclusion");
P("Agriculture must feed a growing population under tightening environmental, climatic and labour constraints, and it must do so sustainably. This chapter has argued that the decisive resource for meeting this challenge is better information, used to make better and more precise decisions. Artificial intelligence, built on the earlier waves of mechanisation, precision agriculture and digitalisation, provides the means to turn the large and growing volume of agricultural data into timely, site-specific and increasingly autonomous action. We have introduced the main families of AI techniques, from expert systems to deep learning and foundation models. We have described a layered, closed-loop architecture for AI-enabled farming and surveyed applications from crop monitoring and disease detection to robotics and supply-chain optimisation. At the same time, the transition is far from complete. Scarce and unrepresentative data, poor generalisation, infrastructure gaps, limited explainability and unresolved questions of data ownership and equity continue to constrain real-world impact, especially for smallholder farmers. Addressing these challenges requires technical innovation together with responsible, inclusive and interdisciplinary practice. The chapters that follow examine each of these themes in depth, with the shared aim of contributing to a smarter, greener and more equitable agricultural future.");

// ---------------------------------------------------------------------------
// Back matter
// ---------------------------------------------------------------------------
const ABBREV = [
  ["AI", "Artificial Intelligence"], ["CNN", "Convolutional Neural Network"], ["DL", "Deep Learning"],
  ["FAO", "Food and Agriculture Organization of the United Nations"], ["GNSS", "Global Navigation Satellite System"],
  ["IoT", "Internet of Things"], ["LoRaWAN", "Long Range Wide Area Network"], ["LSTM", "Long Short-Term Memory"],
  ["ML", "Machine Learning"], ["NB-IoT", "Narrowband Internet of Things"], ["RF", "Random Forest"],
  ["SHAP", "SHapley Additive exPlanations"], ["SVM", "Support Vector Machine"], ["TinyML", "Tiny Machine Learning"],
  ["UAV", "Unmanned Aerial Vehicle"], ["XAI", "Explainable Artificial Intelligence"], ["YOLO", "You Only Look Once"],
];
H1("List of Abbreviations");
for (const [a, f] of ABBREV) {
  body.push(new Paragraph({ spacing: { after: 40 }, tabStops: [{ type: "left", position: 1800 }],
    children: [new TextRun({ text: a, bold: true, font: FONT }), new TextRun({ text: "\t" + f, font: FONT })] }));
}
H1("Acknowledgements");
P("The author thanks the management and the Department of Electronics and Communication Engineering, Velalar College of Engineering and Technology, Erode, for their support and encouragement. [Authors may add funding or other acknowledgements here, or state \"Declared none\".]");
H1("Consent for Publication");
P("Not applicable.");
H1("Conflict of Interest");
P("The author declares no conflict of interest, financial or otherwise.");

H1("References");
order.forEach((key, i) => {
  body.push(new Paragraph({ alignment: AlignmentType.JUSTIFIED, spacing: { after: 80 }, indent: { left: 567, hanging: 567 },
    children: [new TextRun({ text: `[${i + 1}]\t`, font: FONT, size: 22 }), new TextRun({ text: REFS[key], font: FONT, size: 22 })] }));
});
const unused = Object.keys(REFS).filter((k) => !order.includes(k));
if (unused.length) throw new Error("Uncited references: " + unused.join(", "));

// ---------------------------------------------------------------------------
// Front matter (built last so citations above are numbered in body order;
// abstract contains no citations per Bentham guidance)
// ---------------------------------------------------------------------------
const front = [
  new Paragraph({ alignment: AlignmentType.LEFT, spacing: { after: 120 },
    children: [new TextRun({ text: "CHAPTER 1", bold: true, font: FONT, size: 24 })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 120, after: 240 },
    children: [new TextRun({ text: "Introduction to Artificial Intelligence in Agriculture", bold: true, font: FONT, size: 32 })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 120 },
    children: [new TextRun({ text: "M. Nisha Angeline", bold: true, font: FONT, size: 24 }),
      new TextRun({ text: "1,", font: FONT, size: 24, superScript: true }),
      new TextRun({ text: "*", font: FONT, size: 24, superScript: true })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 120 },
    children: [new TextRun({ text: "1", font: FONT, size: 22, superScript: true }),
      new TextRun({ text: " Department of Electronics and Communication Engineering, Velalar College of Engineering and Technology, Thindal, Erode – 638 012, Tamil Nadu, India", italics: true, font: FONT, size: 22 })] }),
  new Paragraph({ alignment: AlignmentType.LEFT, spacing: { after: 360 },
    children: [new TextRun({ text: "*Corresponding author M. Nisha Angeline: ", bold: true, font: FONT, size: 20 }),
      new TextRun({ text: "Department of Electronics and Communication Engineering, Velalar College of Engineering and Technology, Thindal, Erode – 638 012, Tamil Nadu, India; E-mail: nishavlsidesign@gmail.com", font: FONT, size: 20 })] }),
  new Paragraph({ alignment: AlignmentType.JUSTIFIED, spacing: { line: LINE, after: 200 },
    border: { top: { style: BorderStyle.SINGLE, size: 4, color: "808080", space: 6 } },
    children: [new TextRun({ text: "Abstract: ", bold: true, font: FONT }), new TextRun({ font: FONT, text:
      "Agriculture faces the dual challenge of feeding a growing global population while conserving soil, water and biodiversity under a changing climate and growing labour scarcity. Artificial Intelligence (AI) has emerged as a transformative enabler of this transition, converting data from sensors, drones, satellites and machinery into timely, site-specific decisions and autonomous actions. This introductory chapter provides a structured overview of AI in agriculture. It first examines the global drivers that motivate AI adoption and traces the evolution of agricultural technology from mechanisation and precision agriculture to the data-driven paradigms of Agriculture 4.0 and 5.0. The fundamental concepts of AI are then introduced, including knowledge-based systems, machine learning, deep learning, reinforcement learning and emerging foundation models, with a comparative summary of their strengths and limitations. A layered, closed-loop reference architecture for AI-enabled agricultural systems is presented, spanning sensing, connectivity, data management, analytics and actuation. Key application domains are surveyed: crop and soil monitoring, disease and pest detection, precision irrigation and nutrient management, yield prediction, robotics, livestock and supply-chain management. The chapter critically discusses the barriers that limit real-world impact, including data scarcity, poor generalisation, infrastructure and cost constraints, limited explainability and unresolved questions of data ownership and equity. It concludes with emerging research directions and a roadmap of the book." })] }),
  new Paragraph({ alignment: AlignmentType.JUSTIFIED, spacing: { line: LINE, after: 240 },
    border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: "808080", space: 6 } },
    children: [new TextRun({ text: "Keywords: ", bold: true, font: FONT }), new TextRun({ font: FONT, text:
      "Agriculture 4.0, Artificial intelligence, Computer vision, Deep learning, Digital agriculture, Internet of things, Machine learning, Precision agriculture, Smart farming, Sustainable intensification." })] }),
];

const doc = new Document({
  creator: "M. Nisha Angeline",
  title: "Chapter 1: Introduction to Artificial Intelligence in Agriculture",
  styles: {
    default: { document: { run: { font: FONT, size: 24 } } },
    paragraphStyles: [1, 2, 3].map((l) => ({ id: `Heading${l}`, name: `Heading ${l}`, basedOn: "Normal", next: "Normal", quickFormat: true,
      run: { font: FONT, size: 24, bold: true, color: "000000" }, paragraph: { outlineLevel: l - 1, keepNext: true } })),
  },
  numbering: { config: [{ reference: "bullets", levels: [{ level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] }] },
  sections: [{
    properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 20 })] })] }) },
    children: [...front, ...body],
  }],
});

Packer.toBuffer(doc).then((buf) => {
  const out = path.join(__dirname, "Chapter1_Introduction_to_AI_in_Agriculture_Bentham.docx");
  fs.writeFileSync(out, buf);
  console.log("Wrote", out, "with", order.length, "references");
});
