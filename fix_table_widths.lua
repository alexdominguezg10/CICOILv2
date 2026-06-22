-- fix_table_widths.lua
-- Pandoc Lua filter: assign column widths based on column count and header content.
-- Wraps col-1 Code inlines in \texttt{\seqsplit{}} for all table types so that
-- long keys/method names break cleanly instead of overflowing.

local function first_header_text(el)
  local text = ""
  if el.head and el.head.rows and #el.head.rows > 0 then
    local row = el.head.rows[1]
    if row.cells and #row.cells > 0 then
      local cell = row.cells[1]
      for _, block in ipairs(cell.contents) do
        if block.content then
          for _, inline in ipairs(block.content) do
            if inline.text then text = text .. inline.text end
          end
        end
      end
    end
  end
  return text
end

-- Wrap a Code inline in \texttt{\seqsplit{...}} raw LaTeX for breakable monospace.
-- Only escapes underscores; other chars in keys/method names are LaTeX-safe.
local function seqsplit_code(inline)
  if inline.t == "Code" then
    local txt = inline.text:gsub("_", "\\_")
    return pandoc.RawInline("latex", "\\texttt{\\seqsplit{" .. txt .. "}}")
  end
  return inline
end

-- Apply seqsplit_code to all Code inlines in a specific column (0-indexed) of body rows.
local function seqsplit_col(el, col_idx)
  for _, body in ipairs(el.bodies) do
    for _, row in ipairs(body.body) do
      if row.cells and #row.cells > col_idx then
        local cell = row.cells[col_idx + 1]  -- Lua 1-indexed
        local new_contents = {}
        for _, block in ipairs(cell.contents) do
          if block.t == "Para" or block.t == "Plain" then
            local new_inlines = {}
            for _, inline in ipairs(block.content) do
              new_inlines[#new_inlines + 1] = seqsplit_code(inline)
            end
            block.content = new_inlines
          end
          new_contents[#new_contents + 1] = block
        end
        cell.contents = new_contents
      end
    end
  end
  return el
end

local function seqsplit_col1(el)
  return seqsplit_col(el, 0)
end

function Table(el)
  local ncols = #el.colspecs
  local h1 = first_header_text(el)

  if ncols == 2 then
    el.colspecs[1][2] = 0.38
    el.colspecs[2][2] = 0.62
    el = seqsplit_col1(el)

  elseif ncols == 3 then
    if h1 == "Method" or h1 == "Test" then
      -- Method|Signature|Description  or  Test|What|Result
      el.colspecs[1][2] = 0.22
      el.colspecs[2][2] = 0.28
      el.colspecs[3][2] = 0.50
    else
      -- Key|Default|Notes  or  Variable|Units|Default  etc.
      el.colspecs[1][2] = 0.40
      el.colspecs[2][2] = 0.13
      el.colspecs[3][2] = 0.47
    end
    el = seqsplit_col1(el)

  elseif ncols == 4 then
    if h1 == "Variable" then
      -- Variable|Units|Default|Description: defaults are short (just "0"),
      -- so shrink col1+col3 to give col4 more room for long descriptions.
      -- col1 at 35% forces seqsplit breaks on names like fraction_biodegraded_from_water.
      el.colspecs[1][2] = 0.35
      el.colspecs[2][2] = 0.08
      el.colspecs[3][2] = 0.08
      el.colspecs[4][2] = 0.49
    else
      -- Key|Type|Default|Description (§3.1 config table): defaults can be long
      -- plain-text strings like "No_treatment" or "Johansen et al. (2015)".
      el.colspecs[1][2] = 0.40
      el.colspecs[2][2] = 0.08
      el.colspecs[3][2] = 0.20
      el.colspecs[4][2] = 0.32
    end
    el = seqsplit_col1(el)
  end

  return el
end
