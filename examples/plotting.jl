using LaTeXStrings

# Custom plot functions; plots the absolute value of one or many
# wavefunctions with custom axes/legend labels.

# To do:
## Add 2D plot
## Add automatic legend labeling from time-array information.

function my_plot_1D(
  x,
  ψ;
  label,
  kwargs...,
)
  if length(ψ[1])==1
    return plot(
      x,
      abs.(ψ);
      label=label,
      xlabel=L"x",
      ylabel=L"|\psi(x,t)|",
      guidefontsize=16,
      tickfontsize=10,
      legendfontsize=14,
      titlefontsize=16,
      xlims=(-10,10),
      kwargs...,
    )
  else
    return plot(
      x,
      [abs.(i) for i in ψ];
      label=label,
      xlabel=L"x",
      ylabel=L"|\psi(x,t)|",
      guidefontsize=16,
      tickfontsize=10,
      legendfontsize=14,
      titlefontsize=16,
      xlims=(-10,10),
      kwargs...,
    )
  end
end