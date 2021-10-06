class VHDL_Writer :

    def __init__(self, fn):
        self.fn = fn;

    def gen(self, filename):
        with open(filename, "w") as f:
            f.write(f"""
library ieee;
use ieee.std_logic_1164.all;

entity fpr is
    port (
      i_clk :in std_logic;
      i_rst : in std_logic;
      i_seed_data: in std_logic_vector( {len(self.fn) - 1} downto 0);
      output: out std_logic_vector({len(self.fn) - 1} downto 0)
    );
end entity fpr;

architecture run of fpr is

    signal currstate, nextstate:std_logic_vector({len(self.fn) - 1} downto 0);


begin

    statereg: process(i_clk, i_rst)
    begin
        if (i_rst = '1') then
            currstate <= i_seed_data;
        elsif (i_clk = '1' and i_clk'event) then
            currstate <= nextstate;
        end if;
    end process;\n""")
        
            for i in range(len(self.fn) - 1, -1 , -1):
                writestr = ""
                writestr += f"    nextstate({str(i)}) <= "
                for term in self.fn[i]:
                    if type(term) == bool:
                        writestr += f"currstate({str(i)}) XOR "
                    elif len(term) == 1:
                        writestr += f"currstate({str(term[0])}) XOR "
                    else :
                        writestr += "(" + " AND ".join(f"currstate({str(t)})" for t in term) +") XOR "
                writestr = writestr[:-5] + ";\n"
                f.write(writestr);
            f.write("""

    output <= currstate;


end run;



""")
