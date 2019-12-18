require 'nmatrix'
require 'gnuplot'
# 順運動学計算
def f(th)
  px = L1 * Math.cos(th[0]) + L2 * Math.cos(th[0] + th[1])
  py = L1 * Math.sin(th[0]) + L2 * Math.sin(th[0] + th[1])
  return N[[px,py]]
end

# ヤコビ行列の計算
def j(th)
  j00 = -L1 * Math.sin(th[0]) - L2 * Math.sin(th[0] + th[1])
  j01 = -L2 * Math.sin(th[0] + th[1])
  j10 =  L1 * Math.cos(th[0]) + L2 * Math.cos(th[0] + th[1])
  j11 =  L2 * Math.cos(th[0] + th[1])
  return N[[j00,j01],[j10,j11]]
end

# degree から radian への変換
def deg2rad(deg)
  return deg * (Math::PI / 180)
end

# radian から degree への変換
def rad2deg(rad)
  return rad * (180 / Math::PI)
end

L1 = 0.3 # リンク1
L2 = 0.3 # リンク2
k = 0.5 # 収束比例定数
number_of_steps = 20 # 計算回数
r = N[[0.3,0.4]] # 目標位置
th = deg2rad(N[[30,30]]) # 初期関節角度
g = []
orbit_x = []
orbit_y = []
g.push(rad2deg(th))

number_of_steps.times do
  pe = f(th) # 順運動学計算

  # 偏角計算
  dT = NMatrix::BLAS.gemm(j(th).invert,(r-pe),nil,k,0.0,false,:transpose)

  # 関節角度更新
  th = th + dT.transpose

  g.push(rad2deg(th))

  orbit_x.push([0,L1 * Math.cos(th[0]),L1 * Math.cos(th[0]) + L2 * Math.cos(th[0] + th[1])])
  orbit_y.push([0,L1 * Math.sin(th[0]),L1 * Math.sin(th[0]) + L2 * Math.sin(th[0] + th[1])])
end

# 結果の出力
puts rad2deg(th),f(th)

Gnuplot.open do |gp|
  Gnuplot::Plot.new(gp) do |plot|

    plot.output "arm.png"

    x = [0,L1 * Math.cos(th[0]),L1 * Math.cos(th[0]) + L2 * Math.cos(th[0] + th[1])]
    y = [0,L1 * Math.sin(th[0]),L1 * Math.sin(th[0]) + L2 * Math.sin(th[0] + th[1])]

    plot.data << Gnuplot::DataSet.new([x, y]) do |ds|
      ds.with      = "linespoints"
      ds.linewidth = 3
      ds.linecolor = 3
      ds.notitle
    end
  end
end

Gnuplot.open do |gp|
  Gnuplot::Plot.new(gp) do |plot|
    plot.key "under"
    plot.output "theta.png"

    steps = (0..number_of_steps).to_a

    plot.data << Gnuplot::DataSet.new([steps, g.map{|row| row[0]}]) do |ds|
      ds.with      = "linespoints"
      ds.linewidth = 3
      ds.title = "θ1"
    end

    plot.data << Gnuplot::DataSet.new([steps, g.map{|row| row[1]}]) do |ds|
      ds.with      = "linespoints"
      ds.linewidth = 3
      ds.title = "θ2"
    end
  end
end

Gnuplot.open do |gp|
  Gnuplot::Plot.new(gp) do |plot|

    plot.output "arm_orbit.png"

    orbit_x.size.times do |i|
      plot.data << Gnuplot::DataSet.new([orbit_x[i],orbit_y[i]]) do |ds|
        ds.with      = "linespoints"
        ds.linewidth = 1
        ds.notitle
      end
    end
  end
end
