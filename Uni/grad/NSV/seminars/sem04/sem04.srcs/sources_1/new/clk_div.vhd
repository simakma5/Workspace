library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.numeric_std.all;

library UNISIM;
use UNISIM.vcomponents.all;


entity clk_div is
    Port ( CLK100MHZ : in STD_LOGIC;
           CPU_RESETN : in STD_LOGIC;
           clk : out STD_LOGIC);
end clk_div;

architecture Behavioral of clk_div is
    signal cnt: std_logic_vector(31 downto 0);
    signal clk_i: std_logic;
    
    begin
        BUFG_inst: BUFG port map(O => clk, I => clk_i);
        clk_i <= '1' when unsigned(cnt) = 10000000 else '0';
        
        process(CLK100MHZ) begin
            -- or the old fashioned `if CLK100MHZ'event and CLK100MHZ = '1' then`
            if rising_edge(CLK100MHZ) then
                if CPU_RESETN = '0' then
                    cnt <= (others => '0');
                elsif unsigned(cnt) < 10000000 then
                    cnt <= std_logic_vector(unsigned(cnt) + 1);
                else
                    cnt <= (others => '0');
                end if;
            end if;
        end process;

end Behavioral;
