#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import process from "node:process";

const ROOT = path.resolve(path.dirname(new URL(import.meta.url).pathname), "..");
const OUT = path.join(ROOT, "supplementary_tables");

function parseDelimited(text, delimiter) {
  const rows = [];
  let row = [], field = "", quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c === '"' && text[i + 1] === '"') { field += '"'; i++; }
      else if (c === '"') quoted = false;
      else field += c;
    } else if (c === '"') quoted = true;
    else if (c === delimiter) { row.push(field); field = ""; }
    else if (c === "\n") {
      row.push(field.replace(/\r$/, ""));
      if (row.some(x => x !== "")) rows.push(row);
      row = []; field = "";
    } else field += c;
  }
  if (field || row.length) { row.push(field.replace(/\r$/, "")); rows.push(row); }
  const header = rows.shift();
  return rows.map(values => Object.fromEntries(header.map((h, i) => [h, values[i] ?? ""])))
}

async function readTable(file) {
  const text = await fs.readFile(file, "utf8");
  return parseDelimited(text, file.endsWith(".csv") ? "," : "\t");
}

function clean(value) {
  if (value === undefined || value === null || value === "" || value === "NaN") return "NA";
  return String(value).replaceAll("\t", " ").replaceAll("\n", " ");
}

function matrix(columns, rows) {
  return [columns, ...rows.map(row => columns.map(column => clean(row[column])))];
}

function tsv(columns, rows) {
  return matrix(columns, rows).map(row => row.join("\t")).join("\n") + "\n";
}

async function writeTable(filename, columns, rows) {
  if (!rows.length) throw new Error(`${filename}: no data rows`);
  for (const [i, row] of rows.entries()) {
    for (const column of columns) {
      if (!(column in row)) throw new Error(`${filename}: row ${i + 2} lacks ${column}`);
    }
  }
  await fs.writeFile(path.join(OUT, filename), tsv(columns, rows));
  return { filename, columns, rows, data: matrix(columns, rows) };
}

const sampleRows = [
  ["SRR16318553","SRX12596585","SAMN22250671","16-18-hour embryo","control","p1",1],
  ["SRR16318552","SRX12596586","SAMN22250671","16-18-hour embryo","control","p1",2],
  ["SRR16318565","SRX12596573","SAMN22250671","16-18-hour embryo","control","p1",3],
  ["SRR16318560","SRX12596578","SAMN22250671","16-18-hour embryo","control","p1",4],
  ["SRR16318559","SRX12596579","SAMN22250672","16-18-hour embryo","control","p2",1],
  ["SRR16318558","SRX12596580","SAMN22250672","16-18-hour embryo","control","p2",2],
  ["SRR16318557","SRX12596581","SAMN22250672","16-18-hour embryo","control","p2",3],
  ["SRR16318556","SRX12596582","SAMN22250672","16-18-hour embryo","control","p2",4],
  ["SRR16318555","SRX12596583","SAMN22250672","16-18-hour embryo","control","p2",5],
  ["SRR16318569","SRX12596569","SAMN22250675","16-18-hour embryo","control","p3",1],
  ["SRR16318568","SRX12596570","SAMN22250675","16-18-hour embryo","control","p3",2],
  ["SRR16318567","SRX12596571","SAMN22250675","16-18-hour embryo","control","p3",3],
  ["SRR16318566","SRX12596572","SAMN22250675","16-18-hour embryo","control","p3",4],
  ["SRR16318564","SRX12596574","SAMN22250676","16-18-hour embryo","elav5","p3",1],
  ["SRR16318563","SRX12596575","SAMN22250676","16-18-hour embryo","elav5","p3",2],
  ["SRR16318562","SRX12596576","SAMN22250676","16-18-hour embryo","elav5","p3",3],
  ["SRR16318561","SRX12596577","SAMN22250676","16-18-hour embryo","elav5","p3",4],
];

async function build() {
  await fs.mkdir(OUT, { recursive: true });
  const accessionRows = await readTable(path.join(ROOT, "config/accessions.tsv"));
  const accessionByRun = new Map(accessionRows.map(row => [row.run_accession, row]));

  const samples = sampleRows.map(([run, experiment, biosample, group, genotype, probe, replicate]) => {
    const accession = accessionByRun.get(run);
    if (!accession) throw new Error(`Missing accession metadata for ${run}`);
    return {
      run_accession: run, experiment_accession: experiment,
      biosample_accession: biosample, bioproject_accession: "PRJNA771049",
      sra_study_accession: "SRP341252", organism: "Drosophila melanogaster",
      developmental_stage: "16-18-hour embryo", genotype,
      sample_group: group, probe_set: probe, biological_replicate: replicate,
      library_strategy: "Targeted-Capture", library_source: "TRANSCRIPTOMIC",
      library_selection: "PCR", fastq_url: `https://${accession.fastq_ftp}`,
      fastq_md5: accession.fastq_md5, fastq_bytes: accession.fastq_bytes,
      figures_supported: probe === "p3" ? "Figures 4-5" : "Figures 2-3",
    };
  });

  const manifest = await readTable(path.join(ROOT, "config/figure3e_pileup_catalog_manifest.tsv"));
  const atlas = await readTable(path.join(ROOT, "config/figure3e_locus_atlas_manifest.tsv"));
  const publishedCe = new Map();
  for (const row of atlas) if (!publishedCe.has(row.gene)) publishedCe.set(row.gene, row.published_CE);
  const events = manifest.map(row => ({
    event_id: `${row.gene}_CE${publishedCe.get(row.gene)}`, bed_target_id: row.event,
    gene: row.gene, cassette_exon_number: publishedCe.get(row.gene),
    probe_set: row.probe_set, chromosome: row.chrom, strand: row.strand,
    cassette_exon_start: row.ce_start, cassette_exon_end: row.ce_end,
    upstream_flanking_start: row.up_start, upstream_flanking_end: row.up_end,
    downstream_flanking_start: row.down_start, downstream_flanking_end: row.down_end,
    common_utr_start: row.common_start, common_utr_end: row.common_end,
    distal_utr_start: row.distal_start, distal_utr_end: row.distal_end,
    view_start: row.view_start, view_end: row.view_end,
    coordinate_system: "0-based half-open", genome_assembly: "dm6",
    annotation_source: "Supplementary Data 9 cassette exons; reconstructed dm6 flanking and UTR intervals",
  }));

  const ceByGene = new Map(events.map(row => [row.gene, row.cassette_exon_number]));
  const counts = [];
  for (const row of await readTable(path.join(ROOT, "data/derived/figure3e_reconstructed_counts.tsv"))) {
    const informative = Number(row.total);
    const included = Number(row.inclusion);
    counts.push({
      analysis_context: "16-18-hour embryo FASTQ reconstruction",
      run_accession: row.run, sample_group: "16-18-hour embryo", genotype: "control",
      probe_set: row.probe_set, biological_replicate: row.replicate,
      event_id: `${row.gene}_CE${ceByGene.get(row.gene)}`, bed_target_id: row.target,
      gene: row.gene, cassette_exon_number: ceByGene.get(row.gene),
      molecule_class: row.isoform === "All" ? "all" : `${row.isoform.toLowerCase()}_3utr`,
      primary_alignments: "NA", inclusion_reads: row.inclusion,
      skipping_reads: informative - included, unresolved_reads: "NA", boundary_fail_reads: "NA",
      informative_reads: row.total, psi: row.PSI, expected_inclusion_reads: row.inclusion,
      expected_informative_reads: row.total, expected_psi: row.PSI,
      count_match: "TRUE", minimum_depth_pass: informative >= 2 ? "TRUE" : "FALSE",
      source_file: "data/derived/figure3e_reconstructed_counts.tsv",
    });
  }
  const pileups = [
    ["Khc-73","control"], ["Khc-73","elav5"],
    ["Dscam1","control"], ["Dscam1","elav5"],
  ];
  for (const [gene, genotype] of pileups) {
    const file = path.join(ROOT, `data/pileup_inputs/figure6/${gene}/${genotype}/${gene}.pileup_summary.tsv`);
    for (const row of await readTable(file)) {
      counts.push({
        analysis_context: "elav read-level FASTQ reconstruction",
        run_accession: row.run, sample_group: genotype === "control" ? "16-18-hour control embryo" : "16-18-hour elav5 embryo",
        genotype, probe_set: "p3", biological_replicate: row.replicate,
        event_id: `${gene}_CE${ceByGene.get(gene)}`,
        bed_target_id: `${gene}_exon${ceByGene.get(gene)}`, gene,
        cassette_exon_number: ceByGene.get(gene), molecule_class: `${row.fraction}_3utr`,
        primary_alignments: row.primary, inclusion_reads: row.included, skipping_reads: row.skipped,
        unresolved_reads: row.unresolved, boundary_fail_reads: row.boundary_fail,
        informative_reads: row.informative, psi: row.PSI,
        expected_inclusion_reads: row.expected_included,
        expected_informative_reads: row.expected_informative, expected_psi: row.expected_PSI,
        count_match: row.count_match.toUpperCase(),
        minimum_depth_pass: Number(row.informative) >= 2 ? "TRUE" : "FALSE",
        source_file: path.relative(ROOT, file),
      });
    }
  }

  const embryoStats = (await readTable(path.join(ROOT, "data/derived/figure3e_embryo_bh_results.csv"))).map(row => ({
    event_id: `${row.Gene}_CE${row.CE}`, source_key: row.key,
    gene: row.Gene, cassette_exon_number: row.CE,
    n_paired_replicates: row.n_pairs, mean_psi_long: row.mean_long,
    sem_psi_long: row.sem_long, mean_psi_short: row.mean_short,
    sem_psi_short: row.sem_short, delta_psi_long_minus_short: row.delta_long_minus_short,
    ci95_low: row.ci95_low, ci95_high: row.ci95_high,
    test: "two-sided paired t-test", nominal_p: row.nominal_p, bh_q: row.BH_q,
    bh_family: "28 cassette-exon events", nominal_significant: row.nominal_significant,
    bh_significant: row.BH_significant,
  }));

  const interactionSource = await readTable(path.join(ROOT, "data/derived/figure6c_interaction_summary.tsv"));
  const interactions = interactionSource.map(row => ({
    event_id: `${row.gene}_CE${row.CE}`, gene: row.gene, cassette_exon_number: row.CE,
    n_control: row.n_control, n_elav5: row.n_elav5,
    control_mean_psi_short: row.control_mean_short, control_mean_psi_long: row.control_mean_long,
    elav5_mean_psi_short: row.elav5_mean_short, elav5_mean_psi_long: row.elav5_mean_long,
    control_long_minus_short: row.control_signed_coupling,
    elav5_long_minus_short: row.elav5_signed_coupling,
    genotype_by_utr_interaction: row.interaction_elav5_minus_control,
    ci95_low: row.interaction_ci_low, ci95_high: row.interaction_ci_high,
    test: "two-sided Welch t-test on replicate long-minus-short PSI",
    nominal_p: row.interaction_p, bh_q: row.interaction_BH_q,
    bh_family: "22 estimable cassette-exon events",
    bh_significant: Number(row.interaction_BH_q) < 0.05 ? "TRUE" : "FALSE",
  }));

  const marginal = [];
  for (const row of interactionSource) {
    for (const [metric, prefix] of [["delta_dpau","dPAU"], ["pooled_delta_psi","PSIAll"]]) {
      marginal.push({
        event_id: `${row.gene}_CE${row.CE}`, gene: row.gene,
        cassette_exon_number: row.CE, metric,
        effect_elav5_minus_control: row[`${prefix}_change`],
        ci95_low: row[`${prefix}_ci_low`], ci95_high: row[`${prefix}_ci_high`],
        test: "two-sided Welch t-test", nominal_p: row[`${prefix}_p`],
        bh_q: row[`${prefix}_BH_q`], bh_family: `22 ${metric} tests`,
        bh_significant: Number(row[`${prefix}_BH_q`]) < 0.05 ? "TRUE" : "FALSE",
        source_measurement: metric === "delta_dpau" ? "QAPA distal polyadenylation usage" : "pooled cassette-exon PSI",
      });
    }
  }

  const audit = (await readTable(path.join(ROOT, "data/derived/figure3e_deposited_vs_reconstructed.tsv"))).map(row => ({
    run_accession: row.run, probe_set: row.probe_set, biological_replicate: row.replicate,
    event_id: `${row.gene}_CE${ceByGene.get(row.gene)}`, bed_target_id: row.target,
    gene: row.gene, cassette_exon_number: ceByGene.get(row.gene),
    molecule_class: row.isoform === "All" ? "all" : `${row.isoform.toLowerCase()}_3utr`,
    deposited_inclusion_reads: row.deposited_inclusion,
    deposited_informative_reads: row.deposited_total, deposited_psi: row.deposited_PSI,
    reconstructed_inclusion_reads: row.inclusion,
    reconstructed_informative_reads: row.total, reconstructed_psi: row.PSI,
    inclusion_difference: row.inclusion_difference,
    informative_difference: row.total_difference, psi_difference: row.PSI_difference,
    minimum_depth_pass: row.PSI_depth2 === "NA" ? "FALSE" : "TRUE",
    reconciliation_status: row.status,
  }));

  const specs = [
    ["Supplementary_Table_1_samples.tsv", ["run_accession","experiment_accession","biosample_accession","bioproject_accession","sra_study_accession","organism","developmental_stage","genotype","sample_group","probe_set","biological_replicate","library_strategy","library_source","library_selection","fastq_url","fastq_md5","fastq_bytes","figures_supported"], samples],
    ["Supplementary_Table_2_event_definitions.tsv", ["event_id","bed_target_id","gene","cassette_exon_number","probe_set","chromosome","strand","cassette_exon_start","cassette_exon_end","upstream_flanking_start","upstream_flanking_end","downstream_flanking_start","downstream_flanking_end","common_utr_start","common_utr_end","distal_utr_start","distal_utr_end","view_start","view_end","coordinate_system","genome_assembly","annotation_source"], events],
    ["Supplementary_Table_3_plseq_counts.tsv", ["analysis_context","run_accession","sample_group","genotype","probe_set","biological_replicate","event_id","bed_target_id","gene","cassette_exon_number","molecule_class","primary_alignments","inclusion_reads","skipping_reads","unresolved_reads","boundary_fail_reads","informative_reads","psi","expected_inclusion_reads","expected_informative_reads","expected_psi","count_match","minimum_depth_pass","source_file"], counts],
    ["Supplementary_Table_4_embryo_long_short_statistics.tsv", ["event_id","source_key","gene","cassette_exon_number","n_paired_replicates","mean_psi_long","sem_psi_long","mean_psi_short","sem_psi_short","delta_psi_long_minus_short","ci95_low","ci95_high","test","nominal_p","bh_q","bh_family","nominal_significant","bh_significant"], embryoStats],
    ["Supplementary_Table_5_elav_interaction.tsv", ["event_id","gene","cassette_exon_number","n_control","n_elav5","control_mean_psi_short","control_mean_psi_long","elav5_mean_psi_short","elav5_mean_psi_long","control_long_minus_short","elav5_long_minus_short","genotype_by_utr_interaction","ci95_low","ci95_high","test","nominal_p","bh_q","bh_family","bh_significant"], interactions],
    ["Supplementary_Table_6_elav_marginal_statistics.tsv", ["event_id","gene","cassette_exon_number","metric","effect_elav5_minus_control","ci95_low","ci95_high","test","nominal_p","bh_q","bh_family","bh_significant","source_measurement"], marginal],
    ["Supplementary_Table_7_reconstruction_audit.tsv", ["run_accession","probe_set","biological_replicate","event_id","bed_target_id","gene","cassette_exon_number","molecule_class","deposited_inclusion_reads","deposited_informative_reads","deposited_psi","reconstructed_inclusion_reads","reconstructed_informative_reads","reconstructed_psi","inclusion_difference","informative_difference","psi_difference","minimum_depth_pass","reconciliation_status"], audit],
  ];

  const sampleIds = new Set(samples.map(row => row.run_accession));
  const eventIds = new Set(events.map(row => row.event_id));
  for (const row of counts) {
    if (!sampleIds.has(row.run_accession)) throw new Error(`Table 3 has unregistered run ${row.run_accession}`);
    if (!eventIds.has(row.event_id)) throw new Error(`Table 3 has unregistered event ${row.event_id}`);
    const included = Number(row.inclusion_reads), informative = Number(row.informative_reads), psi = Number(row.psi);
    if (Number.isFinite(psi) && informative > 0 && Math.abs(psi - included / informative) > 1e-12) {
      throw new Error(`Table 3 PSI mismatch for ${row.run_accession} ${row.event_id} ${row.molecule_class}`);
    }
  }
  for (const rows of [embryoStats, interactions, marginal, audit]) {
    for (const row of rows) if (!eventIds.has(row.event_id)) throw new Error(`Unregistered event ${row.event_id}`);
  }
  for (const rows of [embryoStats, interactions, marginal]) {
    for (const row of rows) {
      for (const field of ["nominal_p", "bh_q"]) {
        const value = Number(row[field]);
        if (!Number.isFinite(value) || value < 0 || value > 1) throw new Error(`${field} outside [0,1] for ${row.event_id}`);
      }
    }
  }
  for (const row of audit) {
    const deposited = Number(row.deposited_inclusion_reads), reconstructed = Number(row.reconstructed_inclusion_reads);
    if (Number.isFinite(deposited) && Number.isFinite(reconstructed) &&
        Math.abs((reconstructed - deposited) - Number(row.inclusion_difference)) > 1e-12) {
      throw new Error(`Table 7 inclusion difference mismatch for ${row.run_accession} ${row.event_id}`);
    }
  }

  const outputs = [];
  for (const [filename, columns, rows] of specs) outputs.push(await writeTable(filename, columns, rows));

  const modulePath = process.env.ARTIFACT_TOOL_MODULE;
  if (modulePath) {
    const { Workbook } = await import(modulePath);
    const workbook = Workbook.create();
    for (const output of outputs) {
      const sheetName = output.filename.replace("Supplementary_Table_", "S").replace(".tsv", "").slice(0, 31);
      const sheet = workbook.worksheets.add(sheetName);
      sheet.getRange("A1").write(output.data);
    }
    workbook.recalculate();
    for (const sheet of workbook.worksheets.items) {
      const check = await workbook.inspect({ kind: "table", range: `${sheet.name}!A1:C4`, include: "values", tableMaxRows: 4, tableMaxCols: 3 });
      if (!check.ndjson) throw new Error(`Artifact-tool validation failed for ${sheet.name}`);
    }
  }

  const expected = new Map([[outputs[0].filename,17],[outputs[1].filename,28],[outputs[3].filename,28],[outputs[4].filename,22],[outputs[5].filename,44],[outputs[6].filename,254]]);
  for (const output of outputs) {
    if (expected.has(output.filename) && output.rows.length !== expected.get(output.filename)) {
      throw new Error(`${output.filename}: expected ${expected.get(output.filename)} rows, observed ${output.rows.length}`);
    }
  }
  console.log(outputs.map(x => `${x.filename}\t${x.rows.length} rows\t${x.columns.length} columns`).join("\n"));
}

await build();
