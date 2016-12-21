# USAGE:
# ruby double_venn_with_common_genes.rb data/WT0h-WT6h.csv data/KO0h-KO6h.csv results/WT0h-WT6h_vs_KO0h-KO6h_common_genes.xlsx

# ruby double_venn_with_common_genes.rb data/WT0h-WT24h.csv data/KO0h-KO24h.csv results/WT0h-WT24h_vs_KO0h-KO24h_common_genes.xlsx

# ruby double_venn_with_common_genes.rb data/WT6h-WT24h.csv data/KO6h-KO24h.csv results/WT6h-WT24h_vs_KO6h-KO24h_common_genes.xlsx

# Create genes lists with the common genes between 3 DE comparisons.

require 'rubygems'
require 'csv'
require 'axlsx'
require "rinruby" 

ifile1 = ARGV[0]
ifile2 = ARGV[1]
ofile = ARGV[2]

# read the list
counts12_list = Hash.new { |h,k| h[k] = [] }

genes_cond1 = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile1) do |row|
	if row[0] != "id" && row[0] != "" 
		if row[11].to_f >= 1.25 && row[13].to_f <= 0.1
			genes_cond1[row[0]] << 1
		end
	end
end

genes_cond2 = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile2) do |row|
	if row[0] != "id" && row[0] != "" 
		if row[11].to_f >= 1.25 && row[13].to_f <= 0.1
			genes_cond2[row[0]] << 1
		end
	end
end

genes_cond1.each do |gene, count|
	if genes_cond2.has_key?(gene)
		counts12_list[gene] << 1
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

# output
results_xlsx = Axlsx::Package.new
results_wb = results_xlsx.workbook
total_12 = 0

# create list with common peptides between 1,2 conditions
results_wb.add_worksheet(:name => "WT6h-WT24h vs KO6h-KO24h") do |sheet|
	sheet.add_row ["gene", "commons", "WT6h-WT24h count", "KO6h-KO24h count"]
	counts12_list.each do |gene, count|
		sheet.add_row [gene, count.inject(0){|sum,i| sum + i}, genes_cond1[gene].inject(0){|sum,i| sum + i}, genes_cond2[gene].inject(0){|sum,i| sum + i}]
		total_12 += count.inject(0){|sum,i| sum + i}
	end
	sheet.add_row ["totals", total_12, total_1, total_2]
end

# create summary table
results_wb.add_worksheet(:name => "summary") do |sheet|
	sheet.add_row ["total 12", "total 1", "total 2"]
	sheet.add_row [total_12, total_1, total_2]
end

# write xlsx file
results_xlsx.serialize(ofile)

# create venn diagram
R.eval <<EOF
library(VennDiagram)
library(gridExtra)
g = draw.pairwise.venn(area1 = #{total_1}, area2 = #{total_2}, cross.area = #{total_12}, category = c('WT6h-WT24h', 'KO6h-KO24h'), lty = 'blank', fill = c('skyblue', 'mediumorchid'))
pdf('#{ofile}.venn.pdf')
grid.arrange(gTree(children=g), top='DE genes with FC>1.25, qvalue> 0.1')
dev.off()
EOF