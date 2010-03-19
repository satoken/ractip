#!/usr/bin/env ruby

def make_paren_list(x, p_open, p_close)
  res=[]
  xx=x.split('')
  st=[]
  xx.each_index{|i|
    if xx[i]==p_open
      st.push(i)
    elsif xx[i]==p_close
      j=st.pop
      raise "parse error" if j.nil?
      res.push([j, i]) #if j+min_loop<i
    end
  }
  res
end

def get_paren(f)
  open(f) do |fh|
    fh.gets('>')
    p1=fh.gets('>').chomp('>').split(/[\n\r]+/)[2]
    p2=fh.gets('>').chomp('>').split(/[\n\r]+/)[2]
    ex_bp=make_paren_list(p1+p2, '[', ']')
    in_bp1=make_paren_list(p1, '(', ')')
    in_bp2=make_paren_list(p2, '(', ')')
    [ex_bp, in_bp1, in_bp2]
  end
end

def get_acc(tp, ans, res)
  ppv=tp.to_f/res
  sen=tp.to_f/ans
  fval=2*ppv*sen/(ppv+sen)
  [sen, ppv, fval]
end

r=[]
ans=get_paren(ARGV[0])
res=get_paren(ARGV[1])

ex_tp=(ans[0]&res[0]).size
ex_ans=ans[0].size
ex_res=res[0].size
r+=get_acc(ex_tp, ex_ans, ex_res)

in_tp=(ans[1]&res[1]).size+(ans[2]&res[2]).size
in_ans=ans[1].size+ans[2].size
in_res=res[1].size+res[2].size
r+=get_acc(in_tp, in_ans, in_res)

r+=get_acc(ex_tp+in_tp, ex_ans+in_ans, ex_res+in_res)

puts r.map{|v| sprintf("%4.3f", v)}.join(",")
