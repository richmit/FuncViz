#!/usr/bin/env -S ruby
# -*- Mode:ruby; Coding:us-ascii-unix; fill-column:158 -*-
#########################################################################################################################################################.H.S.##
##
# @file      float_diff.rb
# @author    Mitch Richling http://www.mitchr.me/
# @date      2024-08-13
# @brief     File diff comprehending floating point numbers.@EOL
# @std       Ruby 3
# @see       
# @copyright 
#  @parblock
#  Copyright (c) 2024, Mitchell Jay Richling <http://www.mitchr.me/> All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
#
#  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#  3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
#  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  @endparblock
# @filedetails
#
#  Super simple file diff that ignores small differences in floating opts.separator point values.  The regex for floating point values is pretty specific to
#  what my software prints out, so be careful if pressing this script into service for other use cases.
#
#########################################################################################################################################################.H.E.##

require 'optparse'

epsilon   = 1.0e-5
brief     = false
identical = false
opts = OptionParser.new do |opts|
  opts.banner = "Usage: float_diff.rb [options] file1 file2"
  opts.separator ""
  opts.separator "Options:"
  opts.on("-h",         "--help",                   "Show this message")           { puts opts; exit 1;  }
  opts.on("-e epsilon", "--epsilon epsilon",        "Set floating point epsilon")  { |v| epsilon=v.to_f; }
  opts.on("-q",         "--brief",                  "Print just first difference") { |v| brief=true;     }
  opts.on("-s",         "--report-identical-files", "Report identical fiels")      { |v| identical=true; }
  opts.separator ""
  opts.separator "Super simple file diff that ignores small differences in floating"
  opts.separator "point values.  A non-zero exit code is returned if the files are"
  opts.separator "different or an error occurs.  If -a is provided, then all differences"
  opts.separator "are printed.  Otherwise just the fist is printed."
end
opts.parse!(ARGV)

if (ARGV.length < 2) then
  puts("ERROR: Two files must be provided as arguments!")
  exit 2
end

file_size = ARGV.map { |fname|
  begin  
    File::Stat.new(fname).size
  rescue 
    puts("ERROR: Could not stat file argument: '#{fname}'")
    exit 3
  end
}

if (file_size[0] != file_size[1]) then
  puts("#{ARGV[0]} & #{ARGV[1]} have different sizes");
  exit 4
end

file_lines = ARGV.map { |fname|
  open(fname, "r") do |file|
    file.readlines()
  end
}

if (file_lines[0].length != file_lines[1].length) then
  puts("#{ARGV[0]} & #{ARGV[1]} have different line counts");
  exit 5
end

fpre = Regexp.new(/([-+]{0,1}[0-9]\.[0-9]+(e[-+]{0,1}[0-9]+){0,1})/);

num_diffs = 0
file_lines[0].each_index do |idx|
  line_num = idx + 1

  line_floats = [0, 1].map { |i| file_lines[i][idx].scan(fpre).map { |m| m[0].to_f } }

  if (line_floats[0].size != line_floats[1].size) then
    puts("Files have different float counts on line #{line_num}");
    puts("  <<<#{file_lines[0][idx]}")
    puts("  >>>#{file_lines[1][idx]}")
    exit 6 if brief
    num_diffs += 1;
  end

  if (line_floats[0].zip(line_floats[1]).any? { |f0, f1| f0 != f1 }) then
    puts("Files have different float values on line #{line_num}");
    puts("  <<<#{file_lines[0][idx]}")
    puts("  >>>#{file_lines[1][idx]}")
    exit 7 if brief
    num_diffs += 1;
  end

  if (file_lines[0][idx].gsub(fpre, '') != file_lines[1][idx].gsub(fpre, '')) then
    puts("#{ARGV[0]} & #{ARGV[1]} have different non-float content on line #{line_num}");
    puts("  <<<#{file_lines[0][idx]}")
    puts("  >>>#{file_lines[1][idx]}")
    exit 8 if brief
    num_diffs += 1;
  end

end

if (num_diffs > 0) then
  exit 9
else
  puts("#{ARGV[0]} & #{ARGV[1]} are identical!") if identical
  exit 0
end
