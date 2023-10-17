library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity d_ff_tb is
end d_ff_tb;

architecture Behavioral of d_ff_tb is

component d_ff is
    Port ( clk : in STD_LOGIC;
           d : in STD_LOGIC;
           q : out STD_LOGIC;
           rst_n : in STD_LOGIC);
end component;

signal clk: std_logic := '0';
signal d: std_logic := '0';
signal q: std_logic;
signal rst_n: std_logic := '1';

begin

-- instantiate the component d_ff (D flip-flop) and assign its ports
d_ff_i: d_ff
port map(clk => clk, d => d, q => q, rst_n => rst_n);
-- define clock signal behaviour
clk <= not clk after 5ns;

-- define test process
process begin
d <= '0'; rst_n <= '1';
wait for 50ns;
d <= '1'; rst_n <= '1';
wait for 10ns;
d <= '0'; rst_n <= '1';
wait for 10ns;
d <= '1'; rst_n <= '1';
wait for 10ns;
d <= '1'; rst_n <= '1';
wait for 10ns;
d <= '1'; rst_n <= '0';
wait for 10ns;
d <= '1'; rst_n <= '0';
wait for 10ns;
d <= '0'; rst_n <= '0';
wait;
end process;

end Behavioral;
