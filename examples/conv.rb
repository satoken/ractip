#!/usr/bin/env ruby

def conv(t)
  internal=t[2].gsub(/\s/, "").split("")
  external=t[4].gsub(/\s/, "").split("")
  ret=""
  internal.each_index do |i|
    if internal[i]!='.'
      ret+=internal[i]
    elsif external[i]!='.'
      ret+=external[i]
    else
      ret+='.'
    end
  end

  puts ">"+t[0]
  puts t[3].gsub(/(5'-|-3')/, "")
  puts ret
end

l=gets(nil).split(/\n/)
conv(l[0,5])
conv(l[6,5])

