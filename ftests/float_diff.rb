#!/usr/bin/env -S ruby

# Super simple file diff that ignores small differences in floating
# point values.

epsilon = 1.0e-5
allline = true

file_size = ARGV.map { |fname|
  begin  
    File::Stat.new(fname).size
  rescue 
    puts("ERROR: Could not stat file argument: '#{fname}'")
    exit(1)
  end
}

if (file_size[0] != file_size[1]) then
  puts("Files have different sizes");
  exit 3
end

file_lines = ARGV.map { |fname|
  open(fname, "r") do |file|
    file.readlines()
  end
}

if (file_lines[0].length != file_lines[1].length) then
  puts("Files have different line counts");
  exit 4
end

fpre = Regexp.new(/([-+]{0,1}[0-9]\.[0-9]+(e[-+]{0,1}[0-9]+){0,1})/);

file_lines[0].each_index do |idx|
  line_num = idx + 1

  line_floats = [0, 1].map { |i| file_lines[i][idx].scan(fpre).map { |m| m[0].to_f } }

  if (line_floats[0].size != line_floats[1].size) then
    puts("Files have different float counts on line #{line_num}");
    puts("  <<<#{file_lines[0][idx]}")
    puts("  >>>#{file_lines[1][idx]}")
    exit 5 if !allline
  end

  if (line_floats[0].zip(line_floats[1]).any? { |f0, f1| f0 != f1 }) then
    puts("Files have different float values on line #{line_num}");
    puts("  <<<#{file_lines[0][idx]}")
    puts("  >>>#{file_lines[1][idx]}")
    exit 6 if !allline
  end

  if (file_lines[0][idx].gsub(fpre, '') != file_lines[1][idx].gsub(fpre, '')) then
    puts("Files have different non-float content on line #{line_num}");
    puts("  <<<#{file_lines[0][idx]}")
    puts("  >>>#{file_lines[1][idx]}")
    exit 7 if !allline
  end

end
