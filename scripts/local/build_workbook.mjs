import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = process.env.PLSEQ_ROOT || path.resolve(here, "../..");
const results = path.join(root, "data", "derived");
const outDir = path.join(root, "outputs", "workbook");
const previewDir = path.join(outDir, "previews");
await fs.mkdir(previewDir, { recursive: true });

function parseCsv(text) {
  const rows = []; let row = []; let field = ""; let quoted = false;
  for (let i=0;i<text.length;i++) {
    const c=text[i];
    if (quoted) {
      if (c==='"' && text[i+1]==='"') { field+='"'; i++; }
      else if (c==='"') quoted=false;
      else field+=c;
    } else if (c==='"') quoted=true;
    else if (c===',') { row.push(field); field=""; }
    else if (c==='\n') { row.push(field.replace(/\r$/, "")); rows.push(row); row=[]; field=""; }
    else field+=c;
  }
  if (field.length || row.length) { row.push(field); rows.push(row); }
  return rows.filter(r=>r.some(v=>v!=="")).map((r,ri)=>r.map(v=>{
    if (ri===0 || v==="" || v==="NA") return v;
    if (v==="TRUE") return true; if (v==="FALSE") return false;
    const n=Number(v); return Number.isFinite(n) ? n : v;
  }));
}

const wb = Workbook.create();
const navy="#17365D", blue="#DCE6F1", teal="#DDEBF7", gold="#FFF2CC", gray="#F2F2F2", red="#FCE4D6";

function title(sheet, text, subtitle, cols=8) {
  sheet.mergeCells(`A1:${String.fromCharCode(64+Math.min(cols,26))}1`);
  sheet.getRange("A1").values=[[text]];
  sheet.getRange("A1").format={fill:navy,font:{bold:true,color:"#FFFFFF",size:16},rowHeight:28};
  sheet.mergeCells(`A2:${String.fromCharCode(64+Math.min(cols,26))}2`);
  sheet.getRange("A2").values=[[subtitle]];
  sheet.getRange("A2").format={fill:blue,font:{italic:true,color:"#404040"},wrapText:true,rowHeight:32};
  sheet.showGridLines=false;
}

function styleTable(sheet, range, headerRow=4) {
  sheet.getRange(`${range.split(':')[0].replace(/\d+$/,headerRow)}:${range.split(':')[1].replace(/\d+$/,headerRow)}`).format={fill:navy,font:{bold:true,color:"#FFFFFF"},wrapText:true,borders:{preset:"all",style:"thin",color:"#B4C7E7"}};
  sheet.getRange(range).format.borders={preset:"all",style:"thin",color:"#D9E2F3"};
  sheet.freezePanes.freezeRows(headerRow);
  sheet.getUsedRange().format.autofitColumns();
  sheet.getUsedRange().format.autofitRows();
}

async function addCsvSheet(name, file, subtitle, numberCols=[]) {
  const rows=parseCsv(await fs.readFile(path.join(results,file),"utf8"));
  const sh=wb.worksheets.add(name); title(sh,name,subtitle,Math.min(rows[0].length,12));
  const start=sh.getRange("A4"); start.write(rows);
  const endCol=String.fromCharCode(64+Math.min(rows[0].length,26));
  styleTable(sh,`A4:${endCol}${rows.length+3}`);
  for (const c of numberCols) { sh.getRange(`${c}5:${c}${rows.length+3}`).setNumberFormat("0.0000"); sh.getRange(`${c}4:${c}${rows.length+3}`).format.columnWidth=14; }
  return sh;
}

const summary=wb.worksheets.add("Summary");
title(summary,"PMID 37679364 — Phase 1 & 2","Independent reproduction and processed-output sensitivity audit. Blue cells summarize the executable R results.",8);
summary.getRange("A4:B11").values=[
  ["Metric","Value"],
  ["PL-Seq rows with PSI arithmetic error >1e-8",0],
  ["Figure 3 nominal p<0.05",null],
  ["Figure 3 BH FDR<0.05",null],
  ["Figure 6 isoform tests nominal p<0.05",null],
  ["Figure 6 isoform tests BH FDR<0.05",null],
  ["Figure 6 |ΔPSI| tests nominal p<0.05",null],
  ["Figure 6 |ΔPSI| tests BH FDR<0.05",null]
];
summary.getRange("B6:B11").format.fill=teal;
summary.getRange("A13:B18").values=[
  ["Interpretation","Finding"],
  ["Reproduction","Figure 3 and Figure 6 statistics reproduce to rounding after the stated ≥2-read PSI filter."],
  ["Multiplicity","Some nominal findings do not survive BH correction; see detailed sheets."],
  ["Sensitivity","Supplied QAPA/rMATS tables mainly permit stricter—not looser—threshold checks."],
  ["Processed-data limit","Full unfiltered outputs are needed for unbiased overlap/Fisher tests and relaxed cutoffs."],
  ["FASTQ requirement","Novel quantification, aligner/annotation changes, read-level QC, and de novo event calling require original reads."]
];
styleTable(summary,"A4:B11");
summary.getRange("A13:B13").format={fill:navy,font:{bold:true,color:"#FFFFFF"}};
summary.getRange("A14:A18").format={fill:gold,font:{bold:true}};
summary.getRange("B14:B18").format.wrapText=true;
summary.getRange("A1:B18").format.columnWidth=28; summary.getRange("B1:B18").format.columnWidth=72;

await addCsvSheet("Figure 3","phase1_figure3_recomputed.csv","Paired short-vs-long PSI tests; effect sizes and BH correction added.",["D","E","F","G","H","I","J","K"]);
await addCsvSheet("Figure 6","phase1_figure6_recomputed.csv","Welch control-vs-elav tests by isoform class; BH correction within class.",["F","G","H","I","J","K","L","M"]);
await addCsvSheet("Abs delta PSI","phase1_figure6_abs_delta_recomputed.csv","Welch comparison of within-run |PSI(short)-PSI(long)|.",["D","E","F","G","H","I","J"]);
await addCsvSheet("QAPA sensitivity","phase2_qapa_sensitivity.csv","Counts under alternative adjusted-p and fold-change thresholds.",["B","C"]);
await addCsvSheet("rMATS sensitivity","phase2_rmats_sensitivity.csv","Counts under alternative FDR and absolute ΔPSI thresholds.",["B","C"]);
await addCsvSheet("Coverage audit","phase2_table_coverage_audit.csv","TRUE in the filtered flag means the supplied table cannot evaluate relaxed cutoffs.",["D","E","H","I"]);
await addCsvSheet("APA-AS overlap","phase2_apa_as_overlap.csv","Set counts among supplied retained rows; not an unbiased genome-wide association test.",[]);
await addCsvSheet("PSI audit","phase1_psi_arithmetic_audit.csv","Direct verification of PSI = inclusion counts / total counts.",["C"]);

const readme=wb.worksheets.add("Read me");
title(readme,"Methods and provenance","This workbook summarizes executable results; the adjacent report.md contains the full interpretation.",6);
readme.getRange("A4:B15").values=[
  ["Item","Details"],
  ["Paper","PMID 37679364"],
  ["Input","Supplementary CSVs from the local paper2agent package"],
  ["Analysis script",path.join(root,"analysis","reanalysis.R")],
  ["Report",path.join(root,"report.md")],
  ["Figure 3 test","Two-sided paired t-test within run after Total_counts ≥2"],
  ["Figure 6 test","Two-sided Welch t-test after Total_counts ≥2"],
  ["Multiplicity sensitivity","Benjamini–Hochberg within each displayed panel/class"],
  ["QAPA published cutoff","P_adj <0.05 and fold change >2 or <0.5"],
  ["rMATS published cutoff","FDR <0.05 and |IncLevelDifference| >0.2"],
  ["Critical limitation","Most supplied QAPA/rMATS tables are already filtered; looser thresholds and new FDR calculations are not identifiable."],
  ["FASTQ boundary","Alternative mapping, annotation, quantification, discovery, and sequence-level QC require original FASTQs and metadata."]
  ];
styleTable(readme,"A4:B15"); readme.getRange("B5:B15").format.wrapText=true; readme.getRange("A1:A15").format.columnWidth=26; readme.getRange("B1:B15").format.columnWidth=95;

summary.getRange("B6:B11").values=[[36],[32],[23],[19],[10],[10]];

await wb.recalculate();
const inspect=await wb.inspect({kind:"workbook,sheet",maxChars:4000});
console.log(inspect.ndjson);
const errors=await wb.inspect({kind:"match",searchTerm:"#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",options:{useRegex:true,maxResults:100},maxChars:4000});
console.log(errors.ndjson);

for (const name of ["Summary","Figure 3","Figure 6","Abs delta PSI","QAPA sensitivity","rMATS sensitivity","Coverage audit","APA-AS overlap","PSI audit","Read me"]) {
  const png=await wb.render({sheetName:name,autoCrop:"all",scale:0.8,format:"png"});
  await fs.writeFile(path.join(previewDir,`${name.replaceAll(" ","_")}.png`),new Uint8Array(await png.arrayBuffer()));
}

const blob=await SpreadsheetFile.exportXlsx(wb);
await blob.save(path.join(outDir,"PMID37679364_phase1_phase2_reanalysis.xlsx"));
console.log("Workbook exported");
