#!/usr/bin/env -S ruby

epsilon = 1.0e-5

fpre = Regexp.new(/([-+]{0,1}[0-9]\.[0-9]+(e[-+]{0,1}[0-9]+){0,1})/);

begin  
  file0_size = File::Stat.new(ARGV[0]).size
rescue 
  puts("ERROR: Could not stat first file argument: '#{ARGV[0]}'")
  exit(1)
end

begin  
  file1_size = File::Stat.new(ARGV[1]).size
rescue 
  puts("ERROR: Could not stat second file argument: ' #{ARGV[1]}")
  exit(2)
end

if (file0_size != file1_size) then
  puts("Files have diffrent sizes");
  exit 3
end

file0_lines = Array.new
open(ARGV[0], "r") do |file|
  file0_lines = file.readlines()
end

file1_lines = Array.new
open(ARGV[1], "r") do |file|
  file1_lines = file.readlines()
end

if (file0_lines.length != file1_lines.length) then
  puts("Files have diffrent line counts");
  exit 4
end

0.upto(file0_lines.length-1) do |idx|
  line_num = idx + 1

  line0_floats = file0_lines[idx].scan(fpre).map { |m| m[0].to_f }
  line1_floats = file1_lines[idx].scan(fpre).map { |m| m[0].to_f }

  if (line0_floats.size != line1_floats.size) then
    puts("Files have diffrent float counts on line #{line_num}");
    puts("  <<<#{file0_lines[idx]}")
    puts("  >>>#{file1_lines[idx]}")
    exit 5
  end

  line0_floats = file0_lines[idx].scan(fpre).map { |m| m[0].to_f }
  line1_floats = file1_lines[idx].scan(fpre).map { |m| m[0].to_f }

  if (line0_floats.zip(line1_floats).any? { |f0, f1| f0 != f1 }) then
    puts("Files have diffrent float values on line #{line_num}");
    puts("  <<<#{file0_lines[idx]}")
    puts("  >>>#{file1_lines[idx]}")
    exit 6
  end

  if (file0_lines[idx].gsub(fpre, '') != file1_lines[idx].gsub(fpre, '')) then
    puts("Files have diffrent non-float content on line #{line_num}");
    puts("  <<<#{file0_lines[idx]}")
    puts("  >>>#{file1_lines[idx]}")
    exit 7
  end

end
