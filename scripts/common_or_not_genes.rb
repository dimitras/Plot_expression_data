# USAGE:
# ruby common_or_not_genes.rb data/WT0h-WT6h.csv data/KO0h-KO6h.csv WT0-6vsKO0-6

# ruby common_or_not_genes.rb data/WT0h-WT24h.csv data/KO0h-KO24h.csv WT0-24vsKO0-24

# ruby common_or_not_genes.rb data/WT6h-WT24h.csv data/KO6h-KO24h.csv WT6-24vsKO6-24

# Create genes lists with the common genes between 3 DE comparisons.

require 'rubygems'
require 'csv'
require 'axlsx'
require "rinruby" 

ifile1 = ARGV[0]
ifile2 = ARGV[1]
dir = ARGV[2]

condition1 = ifile1.split("/")[1].split(".csv")[0]
condition2 = ifile2.split("/")[1].split(".csv")[0]

# read the lists
genes_cond1 = Hash.new { |h,k| h[k] = [] }
genes_in_1 = Hash.new { |h,k| h[k] = [] }
header_in_1 = nil
CSV.foreach(ifile1) do |row|
	if row[0] == "id"
		header_in_1 = row
	end
	if row[0] != "id" && row[0] != "" 
		if row[11].to_f >= 1.25 && row[13].to_f <= 0.1
			genes_cond1[row[0]] << 1
			genes_in_1[row[0]] << row
		end
	end
end

genes_cond2 = Hash.new { |h,k| h[k] = [] }
genes_in_2 = Hash.new { |h,k| h[k] = [] }
header_in_2 = nil
CSV.foreach(ifile2) do |row|
	if row[0] == "id"
		header_in_2 = row
	end
	if row[0] != "id" && row[0] != "" 
		if row[11].to_f >= 1.25 && row[13].to_f <= 0.1
			genes_cond2[row[0]] << 1
			genes_in_2[row[0]] << row
		end
	end
end

counts12_list = Hash.new { |h,k| h[k] = [] }
genes_cond1.each do |gene, count|
	if genes_cond2.has_key?(gene)
		counts12_list[gene] << 1
	end
end

genes_in_12_from_list1 = Hash.new { |h,k| h[k] = [] }
genes_in_1only = Hash.new { |h,k| h[k] = [] }
genes_in_1.each do |gene, row|
	if genes_in_2.has_key?(gene)
		genes_in_12_from_list1[gene] << row
	else
		genes_in_1only[gene] << row
	end
end

genes_in_12_from_list2 = Hash.new { |h,k| h[k] = [] }
genes_in_2only = Hash.new { |h,k| h[k] = [] }
genes_in_2.each do |gene, row|
	if genes_in_1.has_key?(gene)
		genes_in_12_from_list2[gene] << row
	else
		genes_in_2only[gene] << row
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
CSV.open("results/#{dir}/#{condition1}_only.csv", "wb", {:col_sep => "\t"}) do |csv|
	csv << header_in_1
	genes_in_1only.each do |gene, row|
		csv << row.flatten(2)
	end
end

CSV.open("results/#{dir}/#{condition2}_only.csv", "wb", {:col_sep => "\t"}) do |csv|
	csv << header_in_2
	genes_in_2only.each do |gene, row|
		csv << row.flatten(2)
	end
end

CSV.open("results/#{dir}/#{condition2}in#{condition1}.csv", "wb", {:col_sep => "\t"}) do |csv|
	csv << header_in_1
	genes_in_12_from_list1.each do |gene, row|
		csv << row.flatten(2)
	end
end

CSV.open("results/#{dir}/#{condition1}in#{condition2}.csv", "wb", {:col_sep => "\t"}) do |csv|
	csv << header_in_2
	genes_in_12_from_list2.each do |gene, row|
		csv << row.flatten(2)
	end
end


results_xlsx = Axlsx::Package.new
results_wb = results_xlsx.workbook
total_12 = 0

# create list with common peptides between 1,2 conditions

results_wb.add_worksheet(:name => "commons_counts") do |sheet|
	sheet.add_row ["gene", "commons", "WT counts", "KO counts"]
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
results_xlsx.serialize("results/#{dir}/#{condition1}vs#{condition2}.counts.xlsx")

# create venn diagram
R.eval <<EOF
library(VennDiagram)
library(gridExtra)
g = draw.pairwise.venn(area1 = #{total_1}, area2 = #{total_2}, cross.area = #{total_12}, category = c('#{condition1}', '#{condition2}'), lty = 'blank', fill = c('skyblue', 'mediumorchid'))
pdf('results/#{dir}/#{condition1}vs#{condition2}.venn.pdf')
grid.arrange(gTree(children=g), top='DE genes with FC>1.25, qvalue> 0.1')
dev.off()
EOF