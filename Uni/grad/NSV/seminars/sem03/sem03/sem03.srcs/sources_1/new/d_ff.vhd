library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity d_ff is
    Port ( clk : in STD_LOGIC;
           d : in STD_LOGIC;
           q : out STD_LOGIC;
           rst_n : in STD_LOGIC);
end d_ff;

architecture Behavioral of d_ff is

begin

process(clk) begin
    if rising_edge(clk) then
        if rst_n = '0' then
            q <= '0';
        else
            q <= d;
        end if;
    end if;
end process;

end Behavioral;
