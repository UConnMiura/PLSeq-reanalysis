-- Keep each main figure and its legend together inside one light border in PDF.
-- Other formats retain an ordinary div that can be styled separately.

local function has_class(el, class_name)
  for _, class in ipairs(el.classes) do
    if class == class_name then
      return true
    end
  end
  return false
end

function Div(el)
  if not FORMAT:match("latex") or not has_class(el, "figure-box") then
    return nil
  end

  -- An image with a caption becomes a floating LaTeX figure. Clear the image
  -- caption because the complete figure legend is already a separate paragraph
  -- inside this div and must stay within the same border.
  el = el:walk({
    Image = function(img)
      local path = img.src:gsub("([%%#])", "\\%1")
      return pandoc.RawInline(
        "latex",
        "\\includegraphics[width=\\linewidth,height=0.68\\textheight,keepaspectratio]{" .. path .. "}"
      )
    end
  })

  local content = pandoc.write(pandoc.Pandoc(el.content), "latex")
  local begin_box = [[
\clearpage
\begin{center}
{\color{PaperRule}\fbox{\color{PaperText}\begin{minipage}{0.94\linewidth}
\centering
]]
  local end_box = [[
\end{minipage}}}
\end{center}
\clearpage
]]

  return pandoc.RawBlock("latex", begin_box .. content .. end_box)
end
