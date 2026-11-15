export const GENES=['Glyma.05G103600','Glyma.15G040200','Glyma.20G176100','Glyma.08G163500','Glyma.11G147700'];
export const LABELS={'Glyma.05G103600':'Cortex-1','Glyma.15G040200':'GmHTA12','Glyma.08G163500':'GmAPL-1','Glyma.11G147700':'Epidermis marker'};
export function validateAtlas(d){
 if(!d||d.schema_version!==1||!Array.isArray(d.cells)||!d.cells.length||!d.expression||typeof d.expression!=='object')throw Error('Expected schema_version=1, a nonempty cells array, and an expression object.');
 if(d.cells.length>500000)throw Error('This browser viewer supports up to 500,000 cells. Export treatments separately for larger datasets.');
 const ids=new Set();
 for(const c of d.cells){if(!c||typeof c.id!=='string'||!c.id||ids.has(c.id)||!Number.isFinite(c.x)||!Number.isFinite(c.y)||typeof c.group!=='string'||!c.group||typeof c.celltype!=='string'||!c.celltype)throw Error('Cell IDs must be unique; x/y must be finite; group and celltype must be nonempty strings.');ids.add(c.id);}
 if(!d.cells.some(c=>c.group==='CK'))throw Error('No cells labeled CK were found.');
 if(!Object.keys(d.expression).length)throw Error('No genes were found in expression.');
 for(const [g,v] of Object.entries(d.expression)){if(!Array.isArray(v)||v.length!==d.cells.length||v.some(x=>!Number.isFinite(x)||x<0))throw Error(`${g}: expression length does not match cells, or values are negative/missing.`);}
 return d;
}
export function parseCSV(text){
 text=text.replace(/^\uFEFF/,'');const rows=[];let row=[],field='',quoted=false;
 for(let i=0;i<text.length;i++){const c=text[i];if(c==='"'){if(quoted&&text[i+1]==='"'){field+='"';i++;}else quoted=!quoted;}else if(c===','&&!quoted){row.push(field);field='';}else if((c==='\n'||c==='\r')&&!quoted){if(c==='\r'&&text[i+1]==='\n')i++;row.push(field);if(row.some(x=>x.trim()))rows.push(row);row=[];field='';}else field+=c;}
 if(quoted)throw Error('Unclosed CSV quotation mark.');row.push(field);if(row.some(x=>x.trim()))rows.push(row);
 if(!rows.length)throw Error('CSV is empty.');const header=rows.shift().map(x=>x.trim());const needed=['gene','celltype','treatment','logFC','padj'];if(needed.some(x=>!header.includes(x)))throw Error('Required CSV columns: gene, celltype, treatment, logFC, padj.');
 const seen=new Set();const out=rows.map((r,i)=>{if(r.length!==header.length)throw Error(`CSV row ${i+2} has an incorrect number of columns.`);const o=Object.fromEntries(header.map((h,j)=>[h,r[j].trim()]));if(!o.gene||!o.celltype||!['Heat','Cold','PEG'].includes(o.treatment)||o.logFC===''||o.padj==='')throw Error(`CSV row ${i+2} contains invalid identity fields.`);o.logFC=Number(o.logFC);o.padj=Number(o.padj);if(!Number.isFinite(o.logFC)||!Number.isFinite(o.padj)||o.padj<0||o.padj>1)throw Error(`CSV row ${i+2} has invalid logFC/padj values.`);const key=JSON.stringify([o.gene,o.celltype,o.treatment]);if(seen.has(key))throw Error('Duplicate gene/celltype/treatment rows. Separate time points or models before importing.');seen.add(key);return o;});
 if(!out.length)throw Error('The CSV contains a header but no DE results.');return out;
}
export function sharedCandidates(rows,{fdr=.05,lfc=.5,celltype='all',direction='same'}={}){
 if(!Number.isFinite(fdr)||fdr<0||fdr>1||!Number.isFinite(lfc)||lfc<=0)throw Error('FDR must be between 0 and 1; the absolute log2FC threshold must be positive.');const groups=new Map();
 for(const r of rows){if(celltype!=='all'&&r.celltype!==celltype)continue;const k=JSON.stringify([r.gene,r.celltype]);if(!groups.has(k))groups.set(k,{gene:r.gene,celltype:r.celltype});groups.get(k)[r.treatment]=r;}
 const out=[];for(const g of groups.values()){const v=['Heat','Cold','PEG'].map(t=>g[t]);if(v.some(r=>!r||r.padj>fdr||Math.abs(r.logFC)<lfc))continue;const signs=v.map(r=>Math.sign(r.logFC));const same=signs.every(s=>s===signs[0]);if(direction==='same'&&!same)continue;out.push({...g,maxFDR:Math.max(...v.map(r=>r.padj)),direction:same?(signs[0]>0?'Common up':'Common down'):'Mixed direction'});}return out.sort((a,b)=>a.maxFDR-b.maxFDR||a.gene.localeCompare(b.gene));
}
export function csvCell(v){return '"'+String(v).replaceAll('"','""')+'"';}
export function validateSummaryIndex(d){
 if(!d||d.schema_version!==1||d.type!=='gene_summary_index'||!Array.isArray(d.strata)||!d.strata.length||!d.gene_files||!Object.keys(d.gene_files).length)throw Error('Invalid expression_index.json.');
 for(const s of d.strata)if(typeof s.group!=='string'||typeof s.celltype!=='string'||!Number.isInteger(s.n_cells)||s.n_cells<1)throw Error('Invalid summary strata or cell counts.');
 for(const path of Object.values(d.gene_files))if(typeof path!=='string'||!/^gene_summaries\/[A-Za-z0-9_.-]+\.json$/.test(path))throw Error('Invalid gene summary file path.');
 return d;
}
export function validateGeneSummary(record,index,ncells=null){
 if(!record||!Array.isArray(record.mean)||!Array.isArray(record.pct)||record.mean.length!==index.strata.length||record.pct.length!==index.strata.length)throw Error('Gene summary/strata length mismatch.');
 if(record.mean.some(v=>!Number.isFinite(v)||v<0)||record.pct.some(v=>!Number.isFinite(v)||v<0||v>100))throw Error('Invalid expression mean or detection percentage.');
 if(record.umap){const {i,v}=record.umap;if(!Array.isArray(i)||!Array.isArray(v)||i.length!==v.length||new Set(i).size!==i.length||i.some(k=>!Number.isInteger(k)||k<0||(ncells!==null&&k>=ncells))||v.some(x=>!Number.isFinite(x)||x<=0))throw Error('Invalid sparse single-cell expression vector.');}
 return record;
}
