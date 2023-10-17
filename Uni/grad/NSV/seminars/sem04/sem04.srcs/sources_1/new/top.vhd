library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity top is
    Port ( CLK100MHZ : in STD_LOGIC;
           CPU_RESETN : in STD_LOGIC;
           SW : in STD_LOGIC_VECTOR (1 downto 0);
           LED : out STD_LOGIC_VECTOR (15 downto 0));
end top;

architecture Behavioral of top is
    signal clk: std_logic;
    signal dir_i: std_logic;
    signal en_i: std_logic;
    
    begin
        en_i <= SW(0);
        dir_i <= SW(1);     
        snake_i: entity work.snake port map(
            clk => clk,
            rst_n => CPU_RESETN,
            dir => dir_i,
            en => en_i,
            led => LED);
        clk_div_i: entity work.clk_div port map(
            CLK100MHZ => CLK100MHZ,
            CPU_RESETN => CPU_RESETN,
            clk => clk);

end Behavioral;
