# USAGE:
# ruby triple_venn_with_common_genes.rb data/WT0h-KO0h.csv data/WT6h-KO6h.csv data/WT24h-KO24h.csv results/WT-KOcomparison0.6.24.common_genes_fc1.5pv0.1.xlsx 1.25 0.1

# Create genes lists with the common genes between 3 DE comparisons.

require 'rubygems'
require 'csv'
require 'axlsx'
require "rinruby"

ifile1 = ARGV[0]
ifile2 = ARGV[1]
ifile3 = ARGV[2]
ofile = ARGV[3]
foldchange = ARGV[4].to_f
pvalue = ARGV[4].to_f

# read the list
counts123_list = Hash.new { |h,k| h[k] = [] }
counts12_list = Hash.new { |h,k| h[k] = [] }
counts13_list = Hash.new { |h,k| h[k] = [] }
counts23_list = Hash.new { |h,k| h[k] = [] }
abundance123_list = Hash.new { |h,k| h[k] = [] }
header_in_1 = nil
header_in_2 = nil
header_in_3 = nil

all_genes = Hash.new { |h,k| h[k] = [] }
all_genes1 = Hash.new { |h,k| h[k] = [] }
all_genes2 = Hash.new { |h,k| h[k] = [] }
all_genes3 = Hash.new { |h,k| h[k] = [] }

gene_symbols = {}
fc1 = {}
genes_cond1 = Hash.new { |h,k| h[k] = [] }
genes_cond1_abundance = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile1) do |row|
	if row[0] == "id"
		header_in_1 = row
	end
	if row[0] != "id" && row[0] != "" 
		if row[11].to_f >= foldchange && row[13].to_f <= pvalue
			genes_cond1[row[0]] << 1
			genes_cond1_abundance[row[0]] = row[0..8]
			gene_symbols[row[0]] = row[10]
			fc1[row[0]] = row[11]
			all_genes1[row[0]] = row
		end
	end
end

fc2 = {}
genes_cond2 = Hash.new { |h,k| h[k] = [] }
genes_cond2_abundance = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile2) do |row|
	if row[0] == "id"
		header_in_2 = row
	end
	if row[0] != "id" && row[0] != "" 
		if row[11].to_f >= foldchange && row[13].to_f <= pvalue
			genes_cond2[row[0]] << 1
			genes_cond2_abundance[row[0]] = row[1..8]
			fc2[row[0]] = row[11]
			all_genes2[row[0]] = row
		end
	end
end

fc3 = {}
genes_cond3 = Hash.new { |h,k| h[k] = [] }
genes_cond3_abundance = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile3) do |row|
	if row[0] == "id"
		header_in_3 = row
	end
	if row[0] != "id" && row[0] != "" 
		if row[11].to_f >= foldchange && row[13].to_f <= pvalue
			genes_cond3[row[0]] << 1
			genes_cond3_abundance[row[0]] = row[1..8]
			fc3[row[0]] = row[11]
			all_genes3[row[0]] = row
		end
	end
end

# all_genes1.each do |gene, row|
# 	all_genes[gene] = [row, all_genes2[gene], all_genes3[gene]]
# end

genes_cond1.each do |gene, count|
	if genes_cond2.has_key?(gene) && genes_cond3.has_key?(gene)
		counts123_list[gene] << 1
		abundance123_list[gene] = [genes_cond1_abundance[gene], genes_cond2_abundance[gene], genes_cond3_abundance[gene]]
	end

	if genes_cond2.has_key?(gene)
		counts12_list[gene] << 1
	end

	if genes_cond3.has_key?(gene)
		counts13_list[gene] << 1
	end
end

genes_cond2.each do |gene, count|
	if genes_cond3.has_key?(gene)
		counts23_list[gene] << 1
	end
end

total_1 = 0
genes_cond1.each do |gene, count|
	total_1 += genes_cond1[gene].inject(0){|sum, i| sum + i}
end
total_2 = 0
genes_cond2.each do |gene, count|
	total_2 += genes_cond2[gene].inject(0){|sum, i| sum + i}
end
total_3 = 0
genes_cond3.each do |gene, count|
	total_3 += genes_cond3[gene].inject(0){|sum, i| sum + i}
end

# output
CSV.open("#{ofile}.csv", "wb") do |csv| #, {:col_sep => "\t"}
	csv << header_in_1[0..8].flatten(2) + header_in_2[1..8].flatten(2) + header_in_3[1..8].flatten(2) + ["gene_symbol", "logFC_0h", "logFC_6h", "logFC_24h",]
	abundance123_list.each do |gene, row|
		csv << [row.flatten(2), gene_symbols[gene], fc1[gene], fc2[gene], fc3[gene]].flatten(1)
	end
end

CSV.open("#{ofile}.all_genes.csv", "wb") do |csv|
	csv << header_in_1[0..8].flatten(2) + header_in_2[1..8].flatten(2) + header_in_3[1..8].flatten(2) + ["gene_symbol", "logFC_0h", "logFC_6h", "logFC_24h",]
	all_genes.each do |gene, row|
		csv << [row.flatten(2), gene_symbols[gene], fc1[gene], fc2[gene], fc3[gene]].flatten(1)
	end
end


results_xlsx = Axlsx::Package.new
results_wb = results_xlsx.workbook
total_123 = 0
total_12 = 0
total_13 = 0
total_23 = 0

# create list with common genes between 3 conditions
results_wb.add_worksheet(:name => "common genes in 0,6,24h") do |sheet|
	sheet.add_row ["gene", "commons", "0h count", "6h count", "24h count"]
	counts123_list.each do |gene, count|
		sheet.add_row [gene, count.inject(0){|sum,i| sum + i}, genes_cond1[gene].inject(0){|sum,i| sum + i}, genes_cond2[gene].inject(0){|sum,i| sum + i}, genes_cond3[gene].inject(0){|sum,i| sum + i}]
		total_123 += count.inject(0){|sum,i| sum + i}
	end
	sheet.add_row ["totals", total_123, total_1, total_2, total_3]
end

# create list with common peptides between 1,2 conditions
results_wb.add_worksheet(:name => "common genes in 0,6h") do |sheet|
	sheet.add_row ["gene", "commons", "0h count", "6h count"]
	counts12_list.each do |gene, count|
		sheet.add_row [gene, count.inject(0){|sum,i| sum + i}, genes_cond1[gene].inject(0){|sum,i| sum + i}, genes_cond2[gene].inject(0){|sum,i| sum + i}]
		total_12 += count.inject(0){|sum,i| sum + i}
	end
	sheet.add_row ["totals", total_12, total_1, total_2]
end

# create list with common peptides between 1,3 conditions
results_wb.add_worksheet(:name => "common genes in 0,24h") do |sheet|
	sheet.add_row ["gene", "commons", "0h count", "24h count"]
	counts13_list.each do |gene, count|
		sheet.add_row [gene, count.inject(0){|sum,i| sum + i}, genes_cond1[gene].inject(0){|sum,i| sum + i}, genes_cond3[gene].inject(0){|sum,i| sum + i}]
		total_13 += count.inject(0){|sum,i| sum + i}
	end
	sheet.add_row ["totals", total_13, total_1, total_3]
end

# create list with common peptides between 2,3 conditions
results_wb.add_worksheet(:name => "common genes in 6,24h") do |sheet|
	sheet.add_row ["gene", "commons", "6h count", "24h count"]
	counts23_list.each do |gene, count|
		sheet.add_row [gene, count.inject(0){|sum,i| sum + i}, genes_cond2[gene].inject(0){|sum,i| sum + i}, genes_cond3[gene].inject(0){|sum,i| sum + i}]
		total_23 += count.inject(0){|sum,i| sum + i}
	end
	sheet.add_row ["totals", total_23, total_2, total_3]
end

# create summary table
results_wb.add_worksheet(:name => "summary") do |sheet|
	sheet.add_row ["total 123", "total 12", "total 13", "total 23", "total 1", "total 2", "total 3"]
	sheet.add_row [total_123, total_12, total_13, total_23, total_1, total_2, total_3]
end

# write xlsx file
results_xlsx.serialize(ofile)

# create venn diagram
R.eval <<EOF
library(VennDiagram)
library(gridExtra)
library(extrafont)
# library(ggplot2)
tiff("#{ofile}.venn.tiff", height=800, width=1200, units="px")
g = draw.triple.venn(area1 = #{total_1}, area2 = #{total_2}, area3 = #{total_3}, n12 = #{total_12}, n23 = #{total_23}, n13 = #{total_13}, n123 = #{total_123}, category = c('WT0h-KO0h', 'WT6h-KO6h', 'WT24h-KO24h'), lty = 'blank', fill = c('lightblue', 'blue', 'skyblue'), overrideTriple = TRUE, euler.d = TRUE, scaled=TRUE, cex= rep(3.5, 7), cat.cex = rep(2.5, 3), fontfamily = rep("Arial", 7), cat.fontfamily = rep("Arial", 3))
# grid.arrange(gTree(children=g), top='DE genes with FC>1.25, qvalue> 0.1')
dev.off()

EOF